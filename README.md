Microscopy data needs to first be processed via BacSeg (https://github.com/piedrro/napari-bacseg) to generate segmentations and localise transconjugant foci.
Dependencies:
numpy, pandas, geopandas, pims, scikit-image, scipy, PyQt5, matplotlib

Make sure all the script files are saved in the same directory. 
Set the console working directory to where you are storing these .py files.
To launch the GUI, run Conjugation_GUI.

For each FOV to analyse, the script requires: Donor fluorescence channel images as .tif and the segmentation coordinate outputs from BacSeg in .csv format, these need to be in the same directory. The script also requires the Picasso localisation output file generated via BacSeg. For files to be matched correctly, they need to have the same core naming scheme that consists of the core file name and slice ID. For example "file_name-123_misctext_S1_something.xxx" will have a file name of file_name-123 and a slice id of S1, this will get matched to other files with the same pattern, the Picasso output has a column image_name for each localisation row, these should also follow the abovementioned naming convention. The script also accepts and optional input of the recipient fluorescence channel image as .tif files, these are used to fine-tune the filtering of low-quality localisations for removal of false-positives.

The script saves three dataframes as .csv files. Conjugation_results.csv contains summary info for each FOV analyzed. Conjugation_cells.csv contains the classifiers for each individual cell. Conjugation_localisations.csv contains information about localisations that were assigned as belonging to transconjugants.

The threshold values for assesing localisation quality are currently fixed. 
