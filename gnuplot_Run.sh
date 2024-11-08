#!/bin/bash

# Check if the data file is provided as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <data_file>"
    exit 1
fi

# Define the data file and check if it exists
DATA_FILE="$1"
if [ ! -f "$DATA_FILE" ]; then
    echo "Data file $DATA_FILE does not exist."
    exit 1
fi

# Extract the base name of the data file without extension
BASE_NAME=$(basename "$DATA_FILE" .dat)

# Create a gnuplot script
GNUPLOT_SCRIPT=$(mktemp)
cat <<EOF > $GNUPLOT_SCRIPT
set terminal png
set output 'plots/${BASE_NAME}.png'
set title 'E versus Q^2'
set xlabel 'E_{cm}'
set ylabel 'Q^2'
set palette viridis
plot '$DATA_FILE' using 1:2:3 with points pt 7 ps 1.5 palette notitle
EOF

# Run gnuplot with the script
gnuplot $GNUPLOT_SCRIPT

# Clean up
rm $GNUPLOT_SCRIPT

echo "Plot saved as plots/${BASE_NAME}.png"