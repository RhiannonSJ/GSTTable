# GSTTable
Code to construct a LATEX table using GENIE-generated GST files.

# Running the code

## Inputs

You will need a text file with the format 
    
    <label for table> <path/to/gst/file>

For example

    Nieves2020 /sbnd/data/.../file.txt

With 1 line per input file. If a single file is passed the statistical uncertainty will be written as well.

## Run

The arguments are as follows

    1. Path to the input text file as described above
    2. Detector enumeration of the input file, SBND = 0, MicroBooNE = 1, ICARUS = 2 (for POT scaling)
    3. Path to the desired output file

Run with:

    root -l -b run
    [0] table("textFile.txt",n,"path/to/output/file.tex")

