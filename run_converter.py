import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input file")
args = parser.parse_args()
input_file = args.i

directory = 'converter/'

for filename in os.listdir(directory):
    file_path = os.path.join(directory, filename)
    
    if os.path.isfile(file_path):
        subprocess.run(['python3', file_path, '-i', input_file])
