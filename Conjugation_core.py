import os
import pandas as pd
import numpy as np
import shapely as sp
import geopandas as gpd
import pims
from skimage import draw
from scipy.interpolate import splprep, splev
import re
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed


# Extracts the file name and slice (if applicable)
def extract_file_and_slice(input_string: str) -> tuple[str, str] | tuple[None, None]:
    """
    Extract the file name and slice id from the filepath via regex pattern matching.
    Expected pattern is Filename_123_xxx_S123_xxx
    :param input_string: filepath of segmentation outline .csv file
    :return: file name and slice id as strings
    """
    # define regex pattern to extract file name and slice ID
    pattern = r"(?P<file_id>[\w-]+-\d+)_.*_(?P<slice_id>S\d+)_.*"
    match = re.match(pattern, input_string)
    if match:
        file_id = match.group("file_id")
        slice_id = match.group("slice_id")
        return file_id, slice_id

    else:
        return None, None,


def match_files(file_to_match: str, file_id: str, slice_id: str) -> bool:
    """
    Check if names for segmentation files and localisations match.
    :param file_to_match: filename from the image_name column of the localisations dataframe
    :param file_id: extracted file_id from segmentation file
    :param slice_id: extracted slice_id from segmentation file
    :return: True or False depending on whether the extracted file_id and slice_id match
    """
    pattern = r"(?P<file_id>[\w-]+-\d+)_.*_(?P<slice_id>S\d+)_.*"
    match = re.match(pattern, file_to_match)
    if match and match.group("file_id") == file_id and match.group("slice_id") == slice_id:
        return True
    else:
        return False


def score_poor_locs(locs, quality_parameters: dict) -> int:
    """
    Calculates the quality score of a localisation based on its ellipticity, standard deviation of the PSF, and background.
    :param locs: row from localisations dataframe
    :param quality_parameters: defines the quality thresholds
    :return: quality score of the localisation
    """
    score = 0
    spread_sum = locs.sx + locs.sy
    if locs.ellipticity < quality_parameters['ellipticity_thresholds']['t1']:
        pass
    elif locs.ellipticity < quality_parameters['ellipticity_thresholds']['t2']:
        score -= 1
    elif locs.ellipticity < quality_parameters['ellipticity_thresholds']['t3']:
        score -= 2
    else:
        score -= 5

    if spread_sum < quality_parameters['spread_thresholds']['t1']:
        pass
    elif spread_sum < quality_parameters['spread_thresholds']['t2']:
        score -= 1
    elif spread_sum < quality_parameters['spread_thresholds']['t3']:
        score = -2
    else:
        score -= 5

    if locs.bg < quality_parameters['background_thresholds']['t1']:
        pass
    elif locs.bg < quality_parameters['background_thresholds']['t2']:
        score -= 1
    elif locs.bg < quality_parameters['background_thresholds']['t3']:
        score -= 2
    else:
        score -= 5

    return score


def get_donor_fluorescence(dataframe, fluorescence_image, threshold) -> bool:
    """
    Checks whether a cell is a donor by measuring if its average intracellular fluorescence is above the user defined threshold.
    :param dataframe: dataframe of the cells
    :param fluorescence_image: donor fluorescence channel image
    :param threshold: user defined threshold for what is considered donor fluorescence
    :return:
    """
    mask = np.zeros_like(fluorescence_image, dtype=bool)
    cell_avg_fluorescence_list = []
    x = dataframe.x_new.to_numpy()
    y = dataframe.y_new.to_numpy()
    for i in range(len(dataframe)):
        rr, cc = draw.polygon(y[i], x[i])
        mask[rr, cc] = True
        cell_avg_fluorescence_list.append(np.quantile(fluorescence_image[rr, cc], 0.75))

    bg_fluorescence = np.quantile(fluorescence_image[~mask], 0.10)
    cell_avg_fluorescence_array = np.array(cell_avg_fluorescence_list)
    return (cell_avg_fluorescence_array / bg_fluorescence) > threshold


def recipient_score(first_quartile, bg_fluorescence, threshold):
    cell_ratio = first_quartile / bg_fluorescence
    if cell_ratio < threshold:
        return 1
    elif cell_ratio < threshold * 1.1:
        return 0
    elif cell_ratio < threshold * 1.2:
        return -1
    elif cell_ratio < threshold * 1.3:
        return -2
    elif cell_ratio < threshold * 1.4:
        return -3
    else:
        return -5


def get_recipient_fluorescence(dataframe, fluorescence_image, threshold) -> int:
    """
    Calculates the quality of recipient cell by getting the first quartile of the cellular fluorescence intensity.
    The brighter the cell, the worse it is scored, leading to rejection of localisations within the cell.
    :param dataframe: dataframe of the cells
    :param fluorescence_image: recipient fluorescence channel image
    :param threshold: user defined threshold for scoring recipient quality
    :return: returns the recipient cell quality score which determines whether localisations contained within the cell are considered for analysis
    """
    x = dataframe.x_new.to_numpy()
    y = dataframe.y_new.to_numpy()
    cell_quartile_fluorescence_list = []
    mask = np.zeros_like(fluorescence_image, dtype=bool)
    for i in range(len(dataframe)):
        rr, cc = draw.polygon(y[i], x[i])
        mask[rr, cc] = True
        cell_quartile_fluorescence_list.append(np.quantile(fluorescence_image[rr, cc], 0.25))
    cell_quartile_fluorescence_array = np.array(cell_quartile_fluorescence_list)
    bg_fluorescence = np.quantile(fluorescence_image[~mask], 0.10)
    vfunc = np.vectorize(recipient_score)
    cell_scores = vfunc(cell_quartile_fluorescence_array, bg_fluorescence, threshold)
    print(cell_scores)
    return cell_scores


def classify_cells(row: pd.Series) -> str:
    """
    Assign categorical value to cells based on boolean checks of donor or transconjugant
    :param row: pandas series of cells dataframe
    :return: string categorical cell classification
    """
    if row.is_donor:
        return 'donor'
    elif row.is_transconjugant:
        return 'transconjugant'
    else:
        return 'recipient'


def generate_cells_df(fov_segmentation: pd.DataFrame, fov_localisations: pd.DataFrame,
                      fov_name: str, donor_channel, crop_box,
                      recipient_channel, donor_channel_threshold, recipient_channel_threshold,
                      remove_poor_quality: bool = True):
    # Get cell IDs (~=index) and x y coordinates of the outlines
    work_DF = pd.DataFrame(data={
        'file': fov_name,
        'cell_ID': range(len(fov_segmentation.columns) // 2),
        'x': [np.array(fov_segmentation[f'x[{xi}]']) for xi in range(len(fov_segmentation.columns) // 2)],
        'y': [np.array(fov_segmentation[f'y[{yi}]']) for yi in range(len(fov_segmentation.columns) // 2)]
    })
    # smooth the coordinates via spline interpolation
    residuals = [splprep([work_DF.x.iloc[i][np.logical_not(np.isnan(work_DF.x.iloc[i]))],
                          work_DF.y.iloc[i][np.logical_not(np.isnan(work_DF.y.iloc[i]))]], s=6, k=3) for i in
                 range(work_DF.count().iloc[0])]
    new_points = np.array([splev(residuals[i][1], residuals[i][0]) for i in range(len(residuals))], dtype=object)
    # assign the new x y coordinates of the smoothed outlines
    work_DF['x_new'] = new_points[:, 0]
    work_DF['y_new'] = new_points[:, 1]
    # create geopandas geoseries for storing all the polygon objects
    polygons = gpd.GeoSeries([sp.geometry.Polygon(zip(row.x_new, row.y_new)) for _, row in work_DF.iterrows()])
    # create geopandas geodataframe to store all information about the cells
    cells_gdf = gpd.GeoDataFrame(work_DF[['file', 'cell_ID', 'x_new', 'y_new']], geometry=polygons)
    cells_gdf = cells_gdf[cells_gdf.geometry.within(crop_box)]

    cells_gdf['is_donor'] = get_donor_fluorescence(cells_gdf, donor_channel, donor_channel_threshold)

    # Convert x,y coordinates of localisations to point geometries
    point_locs = gpd.GeoDataFrame(
        fov_localisations, geometry=gpd.points_from_xy(fov_localisations.x, fov_localisations.y)
    )

    if recipient_channel is not None:
        cells_gdf['high_recipient_bg'] = get_recipient_fluorescence(cells_gdf, recipient_channel,
                                                                    recipient_channel_threshold)
        locs_assigned = gpd.sjoin(point_locs, cells_gdf[cells_gdf.is_donor == False])
        locs_assigned = locs_assigned[locs_assigned.quality_score + locs_assigned.high_recipient_bg > -2]
    else:
        locs_assigned = gpd.sjoin(point_locs, cells_gdf[cells_gdf.is_donor == False])
        if remove_poor_quality:
            locs_assigned = locs_assigned[locs_assigned.quality_score > -2]
    cells_gdf['is_transconjugant'] = cells_gdf.cell_ID.isin(locs_assigned.cell_ID)
    cells_gdf['is_recipient'] = ~cells_gdf.is_transconjugant & ~cells_gdf.is_donor
    cells_gdf['cell_type'] = cells_gdf.apply(classify_cells, axis=1)

    return cells_gdf, locs_assigned


def assign_files(donor_tiff_paths: str, localisations_csv_path: str, quality_parameters: dict,
                 cell_fluorescence_thresholds: dict, recipient_tiff_paths: str = None) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_localisations = pd.read_csv(localisations_csv_path)
    fov_name_list = []
    cells_list = []
    locs_list = []
    all_fov_locs_list = []
    recipient_channel = None
    all_recipient_tiff_paths_df = None
    total_recipients = []
    total_transconjugants = []
    total_donors = []
    conjugation_efficiency = []
    transconjugants_per_donor = []
    transconjugants_per_FOV = []
    da_ratio = []
    donor_images = []
    recipient_images = []
    donor_fluorescence_threshold = cell_fluorescence_thresholds['donor_fluorescence_threshold']
    recipient_fluorescence_threshold = cell_fluorescence_thresholds['recipient_fluorescence_threshold']

    if recipient_tiff_paths:
        all_recipient_tiff_paths_df = pd.DataFrame({
            'path': recipient_tiff_paths,
            'basename': [os.path.basename(path) for path in recipient_tiff_paths]})
    else:
        print('No recipient tiff files selected')

    for filepath in donor_tiff_paths:
        segmentation_csv_path = filepath.split('.')[0] + '.csv'
        if os.path.isfile(segmentation_csv_path) and segmentation_csv_path != localisations_csv_path:
            fov_segmentation = pd.read_csv(segmentation_csv_path)
            file_id, slice_id = extract_file_and_slice(os.path.basename(filepath))
            fov_name = f'{file_id}_{slice_id}'
            fov_name_list.append(fov_name)
            fov_localisations = all_localisations[
                all_localisations.image_name.apply(match_files, args=(file_id, slice_id))].copy()
            fov_localisations = fov_localisations[fov_localisations.bg > 0]
            fov_localisations['quality_score'] = fov_localisations.apply(score_poor_locs, args=(quality_parameters,),
                                                                         axis=1)
            fov_localisations['fov_name'] = fov_name
            all_fov_locs_list.append(fov_localisations)
            donor_channel = pims.open(filepath)[0]
            donor_images.append(donor_channel)

            # remove anything near edge of image, in px
            crop_out_width = 10
            x_left, x_right = crop_out_width, donor_channel.shape[1] - crop_out_width - 1
            y_top, y_bottom = crop_out_width, donor_channel.shape[0] - crop_out_width - 1
            crop_box = sp.Polygon([(x_left, y_top), (x_left, y_bottom), (x_right, y_bottom), (x_right, y_top)])

            if all_recipient_tiff_paths_df is not None:
                fov_recipient_tiff_path = all_recipient_tiff_paths_df[
                    all_recipient_tiff_paths_df.basename.apply(match_files, args=(file_id, slice_id,))
                ].path.iloc[0]
                recipient_channel = pims.open(fov_recipient_tiff_path)[0]
                recipient_images.append(recipient_channel)
            else:
                recipient_images.append(np.nan)

            cells_gdf, locs_df = generate_cells_df(fov_segmentation, fov_localisations, fov_name, donor_channel,
                                                   crop_box, recipient_channel, donor_fluorescence_threshold,
                                                   recipient_fluorescence_threshold)
            cells_list.append(cells_gdf)
            locs_list.append(locs_df)

            recipients = (cells_gdf.is_recipient.sum())
            transconjugants = (cells_gdf.is_transconjugant.sum())
            donors = (cells_gdf.is_donor.sum())

            total_recipients.append(recipients)
            total_transconjugants.append(transconjugants)
            total_donors.append(donors)

            conjugation_efficiency.append(get_conventional_conjugation_efficiency(recipients, transconjugants))
            transconjugants_per_donor.append(get_transconjugants_per_donor(transconjugants, donors))
            transconjugants_per_FOV.append(get_transconjugants_per_fov(donors, recipients, transconjugants))
            da_ratio.append(get_cell_ratio(donors, recipients, transconjugants))
            print(fov_name, ' has been processed')

    results_df = pd.DataFrame({
        'file': fov_name_list,
        'total_recipients': total_recipients,
        'total_transconjugants': total_transconjugants,
        'total_donors': total_donors,
        'conjugation_efficiency': conjugation_efficiency,
        'transconjugants_per_FOV': transconjugants_per_FOV,
        'A:D': da_ratio,
        'T:D': transconjugants_per_donor
    })
    images_df = pd.DataFrame({
        'file': fov_name_list,
        'donor_tiff': donor_images,
        'recipient_tiff': recipient_images
    })

    return pd.concat(cells_list), pd.concat(locs_list), results_df, images_df, pd.concat(all_fov_locs_list)


def get_conventional_conjugation_efficiency(total_recipients, total_transconjugants):
    try:
        efficiency = total_transconjugants / (total_recipients + total_transconjugants)
    except ZeroDivisionError:
        return np.nan
    else:
        return efficiency


def get_transconjugants_per_donor(total_transconjugants, total_donors):
    try:
        ratio = total_transconjugants / total_donors
    except ZeroDivisionError:
        return np.nan
    else:
        return ratio


def get_transconjugants_per_fov(total_donors, total_recipients, total_transconjugants):
    try:
        fraction = total_transconjugants / (total_donors + total_recipients + total_transconjugants)
    except ZeroDivisionError:
        return np.nan
    else:
        return fraction


def get_cell_ratio(total_donors, total_recipients, total_transconjugants):
    try:
        ratio = (total_recipients + total_transconjugants) / total_donors
    except ZeroDivisionError:
        return np.nan
    else:
        return ratio

