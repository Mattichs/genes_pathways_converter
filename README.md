# GENE TO PATHWAY CONVERTER FOR MASTRO
## How to run
To run the converter script `run_converter.py`, use the following command:
```
usage: run_converter.py [-h] [-i I]

options:
  -h, --help  show this help message and exit
  -i I        Input file
```
The -i parameter specifies the file to convert.

## Input
The input file must be formatted like a MASTRO input file. For more information, please refer to the [MASTRO documentation](https://github.com/VandinLab/MASTRO).

## Output
Running this script will produce several converted files, depending on the number of converters available in the `converter/` folder. Each converter generates a separate output file, all of which are stored in the `output/` folder. These files can then be used as individual inputs for MASTRO.

## How to add more converter



# Analysis.py
