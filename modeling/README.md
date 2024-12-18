This repository provides MATLAB code for fitting experimental data to different models. The following folders and functions are included:
Folder Structure
1. simple_models

Contains the MATLAB code for fitting experimental data using:

    Activation models.
    Inhibition models.

2. FFL_models

Contains the code for fitting experimental data using the eight feed-forward loop (FFL) models.
3. Model Equations (Models folder)

Each folder includes a Models subfolder where the model equations are written in .def files. These files are required by the Data2Dynamics library for model fitting.
Data Preparation
arrange_data_smad

This script processes the experimental data and creates the Data folder. For each gene:

    A .csv file is generated with the experimental values.
    A .def file is created containing information about the error model.

Fitting Procedure
Setup_smad

This script calls the smadfit function, which performs the fitting process. The results are saved as .mat files.
Library Dependency

The fitting procedure uses the Data2Dynamics library. All input files must be in .def format, as required by the library.
How to Use

    Use arrange_data_smad to prepare the experimental data and create the required Data files.
    Run Setup_smad to fit the models to the data. Results will be saved as .mat files in the appropriate folders.