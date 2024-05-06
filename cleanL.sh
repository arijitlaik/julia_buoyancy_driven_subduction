#!/bin/bash

# Directory where the files and folders are located
DIR="./"  # Current workspace

# Delete the 'markers' directory if it exists
if [ -d "${DIR}/markers" ]; then
    rm -r "${DIR}/markers"
fi

# Delete the .dat files
find $DIR -type f -name "*.dat" -exec rm -f {} \;

# Delete the Proc*.bin files
find $DIR -type f -name "Proc*.bin" -exec rm -f {} \;

# Delete the .vts files
find $DIR -type f -name "*.vts" -exec rm -f {} \;