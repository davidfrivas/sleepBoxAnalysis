# sleepBoxAnalysis

## `compileBoutLengths.m` Documentation

The user has to input a folder containing the `.csv` files (e.g. `/Users/davidrivas/Documents/research/sleep-box/09-18-23`).

```
% Line 12:
mainFolder = '/Users/davidrivas/Documents/research/sleep-box/09-18-23'; % Use current directory, or specify your path
```

It only searches for the files named with `*SB_2sec_*.csv` in the folder, so it doesn't matter if there are also other .csv files in the folder. The files to be analyzed have to be named in this format:
'MM-DD-YYSB_2sec_mouseID.csv' (e.g. `09-18-23SB_2sec_179-6.csv`).

Per mouse (file), it collects all sleep bouts, and bins them by their duration (0s, 4s, 8s, 16s, 32s, 64s, 128s, 256s, 512s). As an example, the 8s bin represents sleep bouts that are greater than or equal to 8 seconds, but less than 16 seconds. The sleep bouts are also grouped in the light phase (6AM-6PM) or dark phase (6PM-6AM).

Then the script calculates statistics across genotypes: the average, std dev, SEM. Then it exports the statistics of individual mice, genotype statistics, and genotype key in separate sheets (`.xls` file). It also outputs some plots in a subfolder `/compiled_plots`.

Sample output:
```
>> compileBoutLengths
Using existing output folder: /Users/davidrivas/Documents/research/sleep-box/09-18-23/compiled_plots
Discovered animal IDs:
  1. 179-6 (from file: 09-18-23SB_2sec_179-6.csv)
  2. 179-7 (from file: 09-18-23SB_2sec_179-7.csv)
  3. 179-8 (from file: 09-18-23SB_2sec_179-8.csv)
  4. 180-3 (from file: 09-18-23SB_2sec_180-3.csv)
  5. 181-4 (from file: 09-18-23SB_2sec_181-4.csv)
  6. 181-5 (from file: 09-18-23SB_2sec_181-5.csv)
  7. 181-6 (from file: 09-18-23SB_2sec_181-6.csv)
  8. 181-7 (from file: 09-18-23SB_2sec_181-7.csv)

Please enter the IDs of wild-type mice, separated by commas:
Wild-type mice: 179-6, 179-8, 180-3, 181-4
Please enter the IDs of mutant mice, separated by commas:
Mutant mice: 179-7, 181-5, 181-6, 181-7
Processing will begin with:
  Wild-type mice: 179-6, 179-8, 180-3, 181-4
  Mutant mice: 179-7, 181-5, 181-6, 181-7

Processing mouse 179-6 (genotype: wild-type)...
Warning: Column headers from the file were modified to make them valid
MATLAB identifiers before creating variable names for the table. The
original column headers are saved in the VariableDescriptions property.
Set 'VariableNamingRule' to 'preserve' to use the original column headers
as table variable names. 
  Found 2998 sleep bouts (1914 light phase, 1084 dark phase)
  Summary for mouse 179-6:
    All sleep bouts: 2998 (binned: [87 467 534 508 508 512 325 56 1])
    Light phase: 1914 (binned: [53 259 286 288 325 380 273 49 1])
    Dark phase: 1084 (binned: [34 208 248 220 183 132 52 7 0])
Processing mouse 179-7 (genotype: mutant)...
...
Creating Excel export with multiple sheets...
Preparing individual mouse data sheet (all sleep)...
  Adding all sleep data for mouse: 179-6
  Adding all sleep data for mouse: 179-8
  Adding all sleep data for mouse: 180-3
  Adding all sleep data for mouse: 181-4
  Adding all sleep data for mouse: 179-7
  Adding all sleep data for mouse: 181-5
  Adding all sleep data for mouse: 181-6
  Adding all sleep data for mouse: 181-7
Preparing individual mouse data sheet (light/dark phases)...
  Adding light/dark phase data for mouse: 179-6
  Adding light/dark phase data for mouse: 179-8
  Adding light/dark phase data for mouse: 180-3
  Adding light/dark phase data for mouse: 181-4
  Adding light/dark phase data for mouse: 179-7
  Adding light/dark phase data for mouse: 181-5
  Adding light/dark phase data for mouse: 181-6
  Adding light/dark phase data for mouse: 181-7
Preparing genotype statistics sheet (all sleep)...
Preparing genotype statistics sheet (light/dark phases)...
Writing data to Excel file: /Users/davidrivas/Documents/research/sleep-box/09-18-23/sleep_bout_analysis.xls
Excel file created successfully.
All figures saved to /Users/davidrivas/Documents/research/sleep-box/09-18-23/compiled_plots
Analysis complete!
```