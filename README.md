# beamformingAnalysis
Flexible analysis tool for beamforming data. Designed to work with TU Delft facilities and MuteSkin database.
 
## Installation
- Clone (or download) this repository
- Add all folders and sub-folders to the MATLAB path

## Usage
- Copy the files in the folder `example` to a folder with beamforming data (either raw .h5 files or .mat files created by the functions in this repository)
- Configure your analysis in `config.xlsx`
- Run the analysis by running `analyse.m` in MATLAB

### Actions
Several analysis actions can be chosen in `config.xlsx`. Below is a typical proceeding:
1. **Inspect raw sample** - Process part of the first file in the folder and check if the beamforming settings give the expected acoustic image. Make some tweaks to the config file, such as refining the integration window
2. **Inspect processed sample** - Reprocess the stored processed sample. This skips the expensive beamforming step, but allows tweaking of the integration window
3. **Batch process raw data** - Process all .h5 files in current folder according to the configurated settings
4. **Convert to database entries** - (For windtunnel measurements of airfoils only) Format the resulting spectra for the MuteSkin database, combining different Reynolds numbers and angles of attack into a single file
