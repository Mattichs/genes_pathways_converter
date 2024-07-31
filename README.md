# GENE TO PATHWAY CONVERTER FOR MASTRO
## Install dependencies
```
pip isntall -r requirements.txt
```

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

## Available conversion types
1. `close_to_root_converter.py`:
- **Description**: This conversion method considers nodes that appear closer to the root node of the tree as more important.
3. `rnd_int_node_converter.py`:
- **Description**: This conversion method considers any internal node as important, regardless of its position relative to the root.

## How to handle more genes and pathways
When multiple trees with unknown genes are used as input for conversion, it is crucial to provide the gene-to-pathway mapping in the `table_gene_path.csv` file located in the `data/` directory. If you wish to abbreviate the names of these pathways, you must also update the `pathways.csv` file, which is also found in the `data/` directory.

## How to add more converter
To add more converter just create a new `update_score()` function with some criterio to update genes score. Then with the new score system, change the `pathwat importance criterio` part. For example, in `close_to_root_converter.py` the criterio check the node with minimum ancestor and if differents nodes have same ancestor, check the which one has maximum descentant.

# Analysis.py
Questo script permette di analizzare i risultati delle varie conversioni rispetto alla rappresentazione originale. E' importante eseguire MASTRO sui dataset convertiti prima di eseguire questa analisi. 

## Input


## Output