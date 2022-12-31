## Getting Started

## Installation
Create a directory to store the repo in by running the following command from your desktop (or wherever you want to store the code)
```
mkdir IRMS
cd IRMS/
```
To clone the repo run the following from the command line in the
```
git clone https://github.com/MDunitz/isotope_ratio_mass_spec.git
```
Create a virtual env
```
python3 -m venv env
```
Activate the virtual env
```
source env/bin/activate
```
Install the required dependencies
```
pip install -r isotope_ratio_mass_spec/requirements.txt
```
## Usage
cd into the directory and activate your virtual env
```
cd IRMS
source env/bin/activate
cd isotope_ratio_mass_spec
```
To label data (assuming the given data is an excel file that meets the expected spec) run the following command
```
python3 -m cli analyze_ms_data PATH_TO_XSLX_DATA_FILE OUTPUT_FILE_NAME
```
To visualize the labled data open the explore_data jupyter notebook under scripts and 

## Citation
