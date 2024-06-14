import matplotlib.colors
import matplotlib.pyplot as plt
from PyQt5.QtGui import QIcon, QPixmap, QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import (
    QApplication, QAction, QMainWindow, QPushButton, QFileDialog,
    QHBoxLayout, QVBoxLayout, QWidget, QMenu, QButtonGroup,
    QRadioButton, QLabel, QCheckBox, QGroupBox, QLineEdit, QComboBox, QFormLayout, QGridLayout)
from PyQt5.QtCore import Qt
import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import sys
import time
import json
from Conjugation_core import *


class MainWindow(QMainWindow):
    def reset_variables(self):
        self.donor_tiff_paths = None
        self.segmentation_csv_paths = None
        self.localisations_csv_path = None
        self.recipient_tiff_paths = None
        self.cells_df = None
        self.locs_df = None
        self.results_df = None
        self.images_df = None
        self.all_fov_locs_df = None

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Conjugation Assay")
        self.setWindowIcon(QIcon(QPixmap("conj_icon.png")))

        # data variables
        self.donor_tiff_paths = None
        self.segmentation_csv_paths = None
        self.localisations_csv_path = None
        self.recipient_tiff_paths = None
        self.cells_df = None
        self.locs_df = None
        self.results_df = None
        self.images_df = None
        self.all_fov_locs_df = None

        menuBar = self.menuBar()
        settingsMenu = menuBar.addMenu("&Settings")
        loc_quality_parameters = QAction("Localisation Quality Parameters", self)
        cell_fluorescence_parameters = QAction("Cell Fluorescence Parameters", self)
        settingsMenu.addAction(loc_quality_parameters)
        settingsMenu.addAction(cell_fluorescence_parameters)
        loc_quality_parameters.triggered.connect(self.open_localisation_parameters_window)
        cell_fluorescence_parameters.triggered.connect(self.open_cell_parameters_window)

        # basic properties of the window

        layout = QVBoxLayout()

        # section for import buttons
        top_layout = QHBoxLayout()
        # import buttons
        self.import_donor_tiff_button = QPushButton("Import Donor Image Channel", self)  # import donor .tiff
        self.import_localisations_csv_button = QPushButton("Import Localisations", self)  # import locs .csv
        self.import_recipient_tiff_button = QPushButton('Import Recipient Image Channel (Optional)', self)
        # connect import buttons to functions
        self.import_donor_tiff_button.clicked.connect(lambda: self.open_multiple_files_dialog('donor'))
        self.import_localisations_csv_button.clicked.connect(self.open_file_dialog)
        self.import_recipient_tiff_button.clicked.connect(lambda: self.open_multiple_files_dialog('recipient'))
        # add import buttons to top_layout
        top_layout.addWidget(self.import_donor_tiff_button)
        top_layout.addWidget(self.import_localisations_csv_button)
        top_layout.addWidget(self.import_recipient_tiff_button)

        # section for processing
        middle_layout = QHBoxLayout()
        self.process_data_button = QPushButton("Process Data", self)
        self.save_results_button = QPushButton("Save Results", self)
        self.save_results_button.setEnabled(False)
        self.plot_results_button = QPushButton("Plot Results", self)
        self.plot_results_button.setEnabled(False)

        # connect widgets to functions
        self.process_data_button.clicked.connect(self.process_data)
        self.save_results_button.clicked.connect(self.save_results)
        self.plot_results_button.clicked.connect(self.plot_results)

        # add widgets to middle layout
        middle_layout.addWidget(self.process_data_button)
        middle_layout.addWidget(self.save_results_button)
        middle_layout.addWidget(self.plot_results_button)

        # last layer for plotting individual FOV results
        bottom_layout = QHBoxLayout()

        self.fov_dropdown = QComboBox(self)
        self.plot_fov_button = QPushButton("Plot FOV", self)
        self.plot_fov_button.setEnabled(False)

        # connect widgets to functions
        self.plot_fov_button.clicked.connect(self.on_plot_fov_button_clicked)
        self.fov_dropdown.currentTextChanged.connect(self.fov_selection_changed)
        bottom_layout.addWidget(self.plot_fov_button)
        bottom_layout.addWidget(self.fov_dropdown)

        layout.addLayout(top_layout)
        layout.addLayout(middle_layout)
        layout.addLayout(bottom_layout)

        self.setCentralWidget(QWidget())
        self.centralWidget().setLayout(layout)
        self.setWindowFlag(Qt.WindowStaysOnTopHint)

    # Methods called by Widgets
    def open_localisation_parameters_window(self):
        self.parameters_window = LocQualityParametersWindow()
        self.parameters_window.show()

    def open_cell_parameters_window(self):
        self.cell_parameters_window = CellFluorescenceWindow()
        self.cell_parameters_window.show()

    def open_multiple_files_dialog(self, button):

        filenames, _ = QFileDialog.getOpenFileNames(
            self, f"Select {button} Channel Files", filter="tif files(*.tif);;All Files (*)"
        )
        if filenames:
            if button == 'donor':
                self.donor_tiff_paths = filenames
                print("Donor files selected")
            else:
                self.recipient_tiff_paths = filenames
        else:
            print("No file selected!")

    def open_file_dialog(self) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select localisations file", filter="csv files (*.csv)")
        if filename:
            self.localisations_csv_path = filename
            print("Localisations file selected")
        else:
            print("No file selected!")

    def process_data(self):
        if self.donor_tiff_paths and self.localisations_csv_path:
            self.fov_dropdown.clear()
            with open('localisation_parameters.json', 'r') as f:
                self.loc_quality_parameters = json.load(f)
            with open('cell_fluorescence_parameters.json', 'r') as f:
                self.cell_fluorescence_parameters = json.load(f)
            start_time = time.time()
            self.cells_df, self.locs_df, self.results_df, self.images_df, self.all_fov_locs_df = assign_files(
                self.donor_tiff_paths,
                self.localisations_csv_path,
                self.loc_quality_parameters,
                self.cell_fluorescence_parameters,
                self.recipient_tiff_paths)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print("Elapsed Time: {:.2f} seconds".format(elapsed_time))
            self.save_results_button.setEnabled(True)
            self.plot_results_button.setEnabled(True)
            self.plot_fov_button.setEnabled(True)
            file_names = self.results_df.file.unique().tolist()
            self.fov_dropdown.addItems(file_names)
        else:
            print("Missing files!")

    def save_results(self):
        save_path = QFileDialog.getExistingDirectory(
            self, caption="Select save directory"
        )
        if save_path:
            self.results_df.to_csv(f'{save_path}/Conjugation_results.csv', index=False)
            self.cells_df.to_csv(f'{save_path}/Conjugation_cells.csv', index=False)
            self.locs_df.to_csv(f'{save_path}/Conjugation_localisations.csv', index=False)
        else:
            print("No save directory selected!")

    def plot_results(self):
        self.plot_results_window = PlotResultsWindow(self.results_df)
        self.plot_results_window.show()

    def on_plot_fov_button_clicked(self):
        self.plot_window = PlotFOVWindow(self.cells_df, self.locs_df, self.images_df, self.selected_fov,
                                         self.all_fov_locs_df)
        self.plot_window.show()

    def fov_selection_changed(self, text):
        self.selected_fov = text


def prepare_histogram(results_df):
    num_FOVs = len(results_df)
    num_graphs = (num_FOVs // 5) + 1 if num_FOVs % 5 != 0 else num_FOVs // 5
    bar_width = 0.2
    FOVs_per_graph = 5
    FOV_numbers = results_df.file
    total_width = num_FOVs * bar_width
    fixed_bar_width = total_width / num_FOVs
    fig, axes = plt.subplots(num_graphs, figsize=(10, 4 * num_graphs))
    axes = np.ravel(axes)  # Flatten the axes array
    fig.suptitle('Cell counts per FOV')
    for i, ax in enumerate(axes):
        start_idx = i * FOVs_per_graph
        end_idx = min((i + 1) * FOVs_per_graph, num_FOVs)

        x = np.arange(start_idx, end_idx)

        # Create bar plots for each type of cell
        ax.bar(x - fixed_bar_width, results_df.total_donors.iloc[start_idx:end_idx], width=fixed_bar_width,
               label='Donors',
               align='center')
        ax.bar(x, results_df.total_recipients.iloc[start_idx:end_idx], width=fixed_bar_width, label='Acceptors',
               align='center')
        ax.bar(x + fixed_bar_width, results_df.total_transconjugants.iloc[start_idx:end_idx], width=fixed_bar_width,
               label='Transconjugants', align='center')

        # Set labels for the x and y axes
        ax.set_ylabel('Counts')

        # Set the x-axis labels to be FOV 1, FOV 2, ..., FOV 8
        ax.set_xticks(x)
        ax.set_xticklabels(FOV_numbers[start_idx:end_idx], fontsize=9)
    fig.legend(labels=['Donors', 'Acceptors', 'Transconjugants'])
    # fig.set_tight_layout(True)
    return fig


def plot_efficiency(results_df):
    fig, ax = plt.subplots()
    fig.suptitle('Conjugation Efficiency')
    ax.set_ylabel("TCJ/(TCJ+Acceptor)")
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.set_xticks([])
    ax.set_ylim(bottom=0)
    ax.boxplot(results_df.conjugation_efficiency)
    fig.set_tight_layout(True)
    return fig


class LocQualityParametersWindow(QWidget):
    def __init__(self):
        super().__init__()

        with open('localisation_parameters.json', 'r') as f:
            self.parameters = json.load(f)

        self.setWindowTitle('Localisation Quality Parameters')
        window_layout = QVBoxLayout()
        self.e_threshold_1 = QLineEdit()
        self.e_threshold_1.setText(str(self.parameters['ellipticity_thresholds']['t1']))
        self.e_threshold_1.setValidator(QDoubleValidator(0.01, 10.00, 2))
        self.e_threshold_2 = QLineEdit()
        self.e_threshold_2.setText(str(self.parameters['ellipticity_thresholds']['t2']))
        self.e_threshold_2.setValidator(QDoubleValidator(0.01, 10.00, 2))
        self.e_threshold_3 = QLineEdit()
        self.e_threshold_3.setText(str(self.parameters['ellipticity_thresholds']['t3']))
        self.e_threshold_3.setValidator(QDoubleValidator(0.01, 10.00, 2))

        self.s_threshold_1 = QLineEdit()
        self.s_threshold_1.setValidator(QDoubleValidator(0.01, 10.00, 2))
        self.s_threshold_1.setText(str(self.parameters['spread_thresholds']['t1']))
        self.s_threshold_2 = QLineEdit()
        self.s_threshold_2.setValidator(QDoubleValidator(0.01, 10.00, 2))
        self.s_threshold_2.setText(str(self.parameters['spread_thresholds']['t2']))
        self.s_threshold_3 = QLineEdit()
        self.s_threshold_3.setText(str(self.parameters['spread_thresholds']['t3']))
        self.s_threshold_3.setValidator(QDoubleValidator(0.01, 10.00, 2))

        self.b_threshold_1 = QLineEdit()
        self.b_threshold_1.setText(str(self.parameters['background_thresholds']['t1']))
        self.b_threshold_1.setValidator(QIntValidator())
        self.b_threshold_2 = QLineEdit()
        self.b_threshold_2.setText(str(self.parameters['background_thresholds']['t2']))
        self.b_threshold_2.setValidator(QIntValidator())
        self.b_threshold_3 = QLineEdit()
        self.b_threshold_3.setText(str(self.parameters['background_thresholds']['t3']))
        self.b_threshold_3.setValidator(QIntValidator())

        ellipticity_threshold_entries = QFormLayout()
        ellipticity_threshold_entries.addRow('Allowed Ellipticity Threshold', self.e_threshold_1)
        ellipticity_threshold_entries.addRow('First Ellipticity Penalty Threshold', self.e_threshold_2)
        ellipticity_threshold_entries.addRow('Second Ellipticity Penalty Threshold', self.e_threshold_3)

        spread_threshold_entries = QFormLayout()
        spread_threshold_entries.addRow('Allowed PSF Spread Threshold', self.s_threshold_1)
        spread_threshold_entries.addRow('First Spread Threshold', self.s_threshold_2)
        spread_threshold_entries.addRow('Second Spread Threshold', self.s_threshold_3)

        background_threshold_entries = QFormLayout()
        background_threshold_entries.addRow('Allowed Localisation Background Threshold', self.b_threshold_1)
        background_threshold_entries.addRow('First Localisation Background Threshold', self.b_threshold_2)
        background_threshold_entries.addRow('Second Localisation Background Threshhold', self.b_threshold_3)

        ellipticity_group_box = QGroupBox('Ellipticity thresholds')
        ellipticity_group_box.setLayout(ellipticity_threshold_entries)
        spread_group_box = QGroupBox('Spread thresholds')
        spread_group_box.setLayout(spread_threshold_entries)
        background_group_box = QGroupBox('Background thresholds')
        background_group_box.setLayout(background_threshold_entries)

        save_button = QPushButton('Save parameters')
        save_button.clicked.connect(self.save_parameters)

        window_layout.addWidget(ellipticity_group_box)
        window_layout.addWidget(spread_group_box)
        window_layout.addWidget(background_group_box)
        window_layout.addWidget(save_button)
        self.setLayout(window_layout)

    def save_parameters(self):
        parameters_dict = {'ellipticity_thresholds': {'t1': float(self.e_threshold_1.text()),
                                                      't2': float(self.e_threshold_2.text()),
                                                      't3': float(self.e_threshold_3.text())},
                           'spread_thresholds': {'t1': float(self.s_threshold_1.text()),
                                                 't2': float(self.s_threshold_2.text()),
                                                 't3': float(self.s_threshold_3.text())},
                           'background_thresholds': {'t1': int(self.b_threshold_1.text()),
                                                     't2': int(self.b_threshold_2.text()),
                                                     't3': int(self.b_threshold_3.text())}
                           }
        with open('localisation_parameters.json', 'w') as outfile:
            json.dump(parameters_dict, outfile)


class CellFluorescenceWindow(QWidget):
    def __init__(self):
        super().__init__()
        with open('cell_fluorescence_parameters.json', 'r') as infile:
            self.localisation_parameters = json.load(infile)

        self.setWindowTitle('Cell Fluorescence Parameters')
        window_layout = QVBoxLayout()
        self.donor_fluorescence_threshold = QLineEdit()
        self.donor_fluorescence_threshold.setText(str(self.localisation_parameters['donor_fluorescence_threshold']))
        self.donor_fluorescence_threshold.setValidator(QDoubleValidator(0.01, 50.00, 2))
        self.recipient_fluorescence_threshold = QLineEdit()
        self.recipient_fluorescence_threshold.setText(
            str(self.localisation_parameters['recipient_fluorescence_threshold']))
        self.recipient_fluorescence_threshold.setValidator(QDoubleValidator(0.01, 50.00, 2))

        cell_fluorescence_thresholds = QFormLayout()
        cell_fluorescence_thresholds.addRow('Donor Fluorescence Threshold', self.donor_fluorescence_threshold)
        cell_fluorescence_thresholds.addRow('Recipient Fluorescence Threshold', self.recipient_fluorescence_threshold)

        save_button = QPushButton('Save')
        save_button.clicked.connect(self.save_fluorescence)

        window_layout.addLayout(cell_fluorescence_thresholds)
        window_layout.addWidget(save_button)

        self.setLayout(window_layout)

    def save_fluorescence(self):
        parameters_dict = {'donor_fluorescence_threshold': float(self.donor_fluorescence_threshold.text()),
                           'recipient_fluorescence_threshold': float(self.recipient_fluorescence_threshold.text())}
        with open('cell_fluorescence_parameters.json', 'w') as outfile:
            json.dump(parameters_dict, outfile)


class PlotResultsWindow(QWidget):
    def __init__(self, results_df):
        super().__init__()
        self.results_df = results_df
        self.setWindowTitle("Conjugation Results")
        global_layout = QHBoxLayout()
        upper_layout = QHBoxLayout()
        upper_canvas = FigureCanvas(prepare_histogram(self.results_df))

        upper_layout.addWidget(NavigationToolbar2QT(upper_canvas, self))
        upper_layout.addWidget(upper_canvas)

        lower_layout = QHBoxLayout()
        lower_canvas = FigureCanvas(plot_efficiency(self.results_df))
        lower_layout.addWidget(NavigationToolbar2QT(lower_canvas, self))
        lower_layout.addWidget(lower_canvas)

        global_layout.addLayout(upper_layout)
        global_layout.addLayout(lower_layout)
        self.setLayout(global_layout)


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # fig.set_tight_layout(True)
        super(MplCanvas, self).__init__(fig)


class PlotFOVWindow(QWidget):
    def __init__(self, cells_df, locs_df, images_df, selected_fov, all_fov_locs):
        super().__init__()
        self.setWindowIcon(QIcon(QPixmap("conj_icon.png")))
        self.setWindowTitle(f"{selected_fov}")
        self.cells_df = cells_df[cells_df.file == selected_fov]
        self.images_df = images_df[images_df.file == selected_fov]
        self.locs_df = locs_df[locs_df.file == selected_fov]
        self.recipient_tiff = self.images_df.recipient_tiff.iloc[0]
        self.all_fov_locs = all_fov_locs[all_fov_locs.fov_name == selected_fov]
        layout = QVBoxLayout()
        tool_layout = QHBoxLayout()

        self.all_locs_checkbox = QCheckBox("Plot all localizations")
        self.all_locs_checkbox.setChecked(True)
        self.all_locs_checkbox.toggled.connect(self.update_plot)
        self.canvas = MplCanvas(self, width=6, height=5, dpi=100)
        self.custom_cmap = matplotlib.colors.ListedColormap(['#1A85FF', '#DDCC77', '#D41159'])
        self.update_plot()

        tool_layout.addWidget(self.all_locs_checkbox)
        tool_layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addLayout(tool_layout)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def update_plot(self):
        self.canvas.axes.cla()
        try:
            self.canvas.axes.imshow(self.recipient_tiff, cmap='Greys_r')
        except TypeError:
            print('No recipient tiff')
        self.cells_df.plot(column='cell_type', legend=True, cmap=self.custom_cmap, facecolor='none', linewidth=2,
                           legend_kwds={'fontsize': 'small', 'loc': 'upper right', 'markerscale': 0.5},
                           ax=self.canvas.axes)
        if self.all_locs_checkbox.isChecked():
            self.canvas.axes.scatter(self.all_fov_locs.x, self.all_fov_locs.y, s=5, color='green', label='all_locs')
        self.canvas.axes.scatter(self.locs_df.x, self.locs_df.y, s=5, color='#D41159')
        self.canvas.draw()


if __name__ == '__main__':
    plt.ioff()
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())



