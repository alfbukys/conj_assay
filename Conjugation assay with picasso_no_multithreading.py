import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shapely as sp
from scipy.io import loadmat
import geopandas as gpd
import pims
from skimage import draw
from scipy.interpolate import splprep, splev
import seaborn as sns
import tkinter as tk
from tkinter import filedialog
import pickle
import re
import time

    
# Extracts the file name and slice (if applicable)
def extract_file_and_slice(input_string):
    # define regex pattern to extract file name and slice ID
    pattern = r"(?P<file_id>[\w-]+-\d+)_.*_(?P<slice_id>S\d+)_.*"
    match = re.match(pattern, input_string)
    if match:
        file_id = match.group("file_id")
        slice_id = match.group("slice_id")
        return file_id, slice_id

    else:
        return None, None, 

    
def match_files(file_to_match, file_id, slice_id,):
    pattern = r"(?P<file_id>[\w-]+-\d+)_.*_(?P<slice_id>S\d+)_.*"
    match = re.match(pattern, file_to_match)
    if match and match.group("file_id") == file_id and match.group("slice_id") == slice_id:
        return True
    else:
        return False
    
    
def file_importer (useIndex=False):
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    
    # Raise the main window to the top (optional)
    root.attributes('-topmost', True)
    root.update()
    
    # Open a directory dialog to select a directory
    directory = filedialog.askdirectory(title="Select a Directory")
    picasso_locs_path = filedialog.askopenfilename(title='Select localisation file', filetypes=[('csv files', '*.csv')])
    all_locs_df = pd.read_csv(picasso_locs_path)
    # Close the hidden main window (optional)
    root.destroy() 
    
    if directory and picasso_locs_path:
        file_dict = {} #dictionary to store imported data with a key for the donor channel file name
        # output_dict = {} #dictionary to store processed data with a key for the donor channel file name
        index = 0
        for filename in os.listdir(directory):
            if filename.endswith(".tif") or filename.endswith(".tiff"):
                base_filename = os.path.splitext(filename)[0]
                file_id, slice_id = extract_file_and_slice(base_filename)
                csv_filename = base_filename + ".csv"
                csv_path = os.path.join(directory, csv_filename)
                if os.path.isfile(csv_path) and csv_path != picasso_locs_path:
                    tiff_path = os.path.join(directory, filename)
                    if useIndex:
                        all_fov_locs = all_locs_df[all_locs_df.frame == index]
                    else:
                        all_fov_locs = all_locs_df[all_locs_df.image_name.apply(match_files, args =(file_id, slice_id))]
                        # all_fov_locs = all_locs_df[all_locs_df.image_name.str.contains(file_id) &
                        #                             all_locs_df.image_name.str.contains(slice_id)]
                    simple_name = f'{file_id}_{slice_id}'
                    input_data = Data_input(tiff_path, csv_path, all_fov_locs, simple_name)
                    file_dict[simple_name] = input_data
                    # output_dict[base_filename] = classify_cells(input_data.segmentation, all_fov_locs, input_data.image, base_filename)
                    print(f'{simple_name} has been imported')
                    index +=1
    return (file_dict, directory)
    

def generate_outline (segmentation, isOufti):
    if isOufti == True:
        temp_handling =  [[item['model'].flat[0] for item in cell] for cell in segmentation['cellList'][0][0][0]]
        temp_handling = list(filter(None, temp_handling))
        Polygons = gpd.GeoDataFrame({'cell_ID' : [i for i in range(len(temp_handling))],'geometry' : [sp.geometry.Polygon(item[0]) for item in temp_handling] })
    else:
        #work_DF prepares csv data for geopandas and shapely
        work_DF = pd.DataFrame(data={'cell_ID' : [i for i in range(len(segmentation.columns)//2)],
                                     'x' : [np.array(segmentation['x[{xi}]'.format(xi=i)]) for i in range(len(segmentation.columns)//2)],
                                     'y' : [np.array(segmentation['y[{yi}]'.format(yi=i)]) for i in range(len(segmentation.columns)//2)]})
        #ressies stores the cubic interpolation of cell coordinates
        ressies = [splprep([work_DF.x.iloc[i][np.logical_not(np.isnan(work_DF.x.iloc[i]))], work_DF.y.iloc[i][np.logical_not(np.isnan(work_DF.y.iloc[i]))]], s=6, k=3) for i in range(work_DF.count().iloc[0])]
        #new points generates new coordinates from the interpolation
        new_points = np.array([splev(ressies[i][1], ressies[i][0]) for i in range(len(ressies))], dtype=object)

        #replace old points with the new ones in the dataframe
        work_DF['x_new'] = new_points[:,0]
        work_DF['y_new'] = new_points[:,1]

        #generate polygons geoDataFrame
        Polygons = gpd.GeoDataFrame({'cell_ID' : work_DF.cell_ID, 
                                      'geometry' : [sp.geometry.Polygon(zip(row.x_new, row.y_new)) for index, row in work_DF.iterrows()]
                                      })
    return(Polygons)


def get_fluorescence(Outline, fluorescence_image):
    mask = np.zeros_like(fluorescence_image, dtype=bool)
    rr, cc = draw.polygon(Outline.geometry.iloc[0].exterior.xy[1], Outline.geometry.iloc[0].exterior.xy[0])
    mask[rr, cc] = True
    masked_image = np.copy(fluorescence_image)
    masked_image[~mask] = 0
    avg_fluor = np.mean(masked_image[mask])
    return(avg_fluor)

#removes dead pixel localisations and other dubious locs, change max_blur if you find it removes too many true localisations
def clean_locs(locs):
    max_blur = 20
    cleaned_locs = locs[ ( 0.8 < locs.sx + locs.sy) & 
                        (locs.bg < 1800) & 
                        (locs.net_gradient >= 4500) &
                        (locs.sx + locs.sy < max_blur) &
                        (locs.bg > 0)
                        ]
    return(cleaned_locs)
    

def classify_cells(input_object):
    segmentation = input_object.segmentation
    all_locs = clean_locs(input_object.all_locs)
    donor_fluorescence_image = input_object.image
    base_filename = input_object.base_filename
    isOufti = input_object.base_filename
    Outlines = generate_outline(segmentation, isOufti)#get all outlines of all cells
    Point_locs = gpd.GeoDataFrame(
        all_locs, geometry = gpd.points_from_xy(all_locs.x, all_locs.y)
        )
    total_donors, total_recipients, total_transconjugants = 0, 0, 0 #initiate counting of cell types
    Cells_dict = {} #dictionary to store info about individual cells
    crop_out_width = 10
    x_left, x_right = crop_out_width, 2752 - crop_out_width - 1
    
    y_top, y_bottom = crop_out_width, 2208 - crop_out_width - 1
    crop_box = sp.Polygon([(x_left, y_top), (x_left, y_bottom), (x_right, y_bottom), (x_right, y_top)])
    for i in range(len(Outlines)): #we will be iterating through each cell individually
        cell_outline = Outlines.loc[Outlines.cell_ID == i] #get outline of the cell
        if not cell_outline.values[0, 1].within(crop_box) :
            continue
        fluor_signal = get_fluorescence(cell_outline, donor_fluorescence_image) #quantify fluorescence signal within the cell
        if fluor_signal > 500: 
            cell_type = 'donor' #if signal above threshold, cell is donor
            total_donors += 1
        elif (Point_locs.within(cell_outline.geometry.iloc[0], align=False)).any(): 
            cell_type = 'transconjugant' #otherwise check if contain foci
            total_transconjugants += 1
        else: 
            cell_type = 'recipient' #if empty for both, cell is a recipient
            total_recipients += 1
        Cells_dict[i] = Cell_info(cell_outline, cell_type, Point_locs.loc[Point_locs.within(cell_outline.geometry.iloc[0])])
    conventional_conjugation_efficiency = get_conventional_conjugation_efficiency(total_recipients, total_transconjugants)
    transconjugants_per_FOV = get_transconjugants_per_fov(total_donors, total_recipients, total_transconjugants)
    donor_to_recipient_ratio = get_cell_ratio(total_donors, total_recipients, total_transconjugants)
    cell_counts = pd.DataFrame(data={'file':[base_filename],
                                     'donors':[total_donors], 
                                     'recipients':[total_recipients], 
                                     'transconjugants':[total_transconjugants],
                                     'conventional_conjugation_efficiency':[conventional_conjugation_efficiency],
                                     'transconjugants_per_FOV':[transconjugants_per_FOV],
                                     'donor_to_recipient_ratio':[donor_to_recipient_ratio]})
    return Data_output(Cells_dict, cell_counts)


def process_data(key, data):
    result = classify_cells(data)    
    print(f'{key} has been processed')
    return result

    
def get_conventional_conjugation_efficiency(total_recipients, total_transconjugants):
    efficiency = total_transconjugants/(total_recipients + total_transconjugants)
    return efficiency

def get_transconjugants_per_fov(total_donors, total_recipients, total_transconjugants):
    fraction = total_transconjugants / (total_donors + total_recipients + total_transconjugants)
    return fraction

def get_cell_ratio(total_donors, total_recipients, total_transconjugants):
    ratio = (total_recipients + total_transconjugants)/total_donors
    return(ratio)

def generate_summary(Data_output):
    results_data = {
        key : value.results for key, value in Data_output.items()
        }
    results_df = pd.concat([values for values in results_data.values()])
    avg_cells = results_df.mean(numeric_only=True).to_frame().T
    return(results_df, avg_cells)

def plot_cell_counts_by_FOV(results_df, directory):
    # Number of fields of view
    num_FOVs = len(results_df)

    # Determine how many graphs are needed (e.g., with a maximum of 8 bars per graph)
    num_graphs = (num_FOVs // 5) + 1 if num_FOVs % 5 != 0 else num_FOVs // 5

    # Create an array of FOV numbers (e.g., FOV 1, FOV 2, ..., FOV 8)
    FOV_numbers = results_df.file

    # Set the width of each bar
    bar_width = 0.2

    # Set the number of FOVs per graph
    FOVs_per_graph = 5

    # Calculate the total width of all bars
    total_width = num_FOVs * bar_width

    # Set the width for each individual bar
    fixed_bar_width = total_width / num_FOVs

    # Create a figure with subplots
    fig, axes = plt.subplots(num_graphs, figsize=(10, 4 * num_graphs))
    axes = np.ravel(axes)  # Flatten the axes array
    fig.suptitle('Cell counts per FOV')
    for i, ax in enumerate(axes):
        start_idx = i * FOVs_per_graph
        end_idx = min((i + 1) * FOVs_per_graph, num_FOVs)

        x = np.arange(start_idx, end_idx)

        # Create bar plots for each type of cell
        ax.bar(x - fixed_bar_width, results_df.donors.iloc[start_idx:end_idx], width=fixed_bar_width, label='Donors', align='center')
        ax.bar(x, results_df.recipients.iloc[start_idx:end_idx], width=fixed_bar_width, label='Acceptors', align='center')
        ax.bar(x + fixed_bar_width, results_df.transconjugants.iloc[start_idx:end_idx], width=fixed_bar_width, label='Transconjugants', align='center')

        # Set labels for the x and y axes
        ax.set_xlabel('Fields of View')
        ax.set_ylabel('Counts')


        # Set the x-axis labels to be FOV 1, FOV 2, ..., FOV 8
        ax.set_xticks(x)
        ax.set_xticklabels(FOV_numbers[start_idx:end_idx], fontsize=9)
        ax.legend()
    fig.savefig(f'{directory}/cell_counts_by_FOV.png', dpi=300)

def plot_conjugation_efficiencies(results_df, directory):
    # Create a subplot
    fig, ax = plt.subplots()
    sns.violinplot(y=results_df.conventional_conjugation_efficiency, 
                   color='#D81B60', ax=ax).set(title='Conventional conjugation efficiency')
    # Set labels for the plot
    ax.set_ylabel("TCJ/(TCJ+Acceptor)")
    fig.savefig(f'{directory}/Conventional_conjugation_efficiency.png', dpi=300)
    # ax.set_ylim(0)
    
    # Set labels for the plot
    ax.set_ylabel("TCJ/Donor")
    fig.savefig(f'{directory}/Transconjugants_per_donor.png', dpi=300)
    # ax.set_ylim(0)
    
    # Set labels for the plot
    ax.set_ylabel("TCJ/(Donor x (Recipient+TCJ)")
    fig.savefig(f'{directory}/Normalized_conjugation_efficiency.png', dpi=300)
    # ax.set_ylim(0)
    

class Data_input:
    def __init__(self, donor_path, segmentation_path, all_locs, base_filename, isOufti=False):
        self.isOufti = isOufti
        if isOufti:
            self.segmentation = loadmat(segmentation_path)
        else:
            self.segmentation = pd.read_csv(segmentation_path)
        self.image = pims.open(donor_path)[0]
        self.all_locs = all_locs
        self.base_filename = base_filename
        
class Cell_info: #stores values about individual cells
    def __init__(self, outline, cell_type, locs):
        self.outline = outline #stores the geopandas dataframe with Shapely polygons
        self.cell_type = cell_type #determines the cell type (donor, recipient or transconjugant)
        self.locs = locs #stores the localisations
        
class Data_output:
    def __init__(self, cell_data, results):
        self.cell_data = cell_data
        self.results = results
        


Input_dict, directory = file_importer() #set useIndex to true if regex matching breaks, won't work if concatinating from multiple sessions
start_time = time.time()
Output_dict = {key : process_data(key, data) for (key, data) in Input_dict.items()}
end_time = time.time()
elapsed_time = end_time - start_time
print("Elapsed Time: {:.2f} seconds".format(elapsed_time))   

results, summary_df = generate_summary(Output_dict)
results_df = results.sort_values(by='file')



 
plot_cell_counts_by_FOV(results_df, directory)
plot_conjugation_efficiencies(results_df, directory)

results_df.to_csv(f'{directory}/Results.csv', index=False)
summary_df.to_csv(f'{directory}/Summary.csv', index=False)
with open(f'{directory}/Processed_data.pkl', 'wb') as f:
    pickle.dump((Output_dict, results_df, summary_df), f)