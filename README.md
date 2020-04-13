# Phast User Event

This repository stores all the code to use the Phast User Event producing ROOT TTrees for Multiplicity analysis.

## User Event

The [UserEvent20](https://github.com/nipierre/Phast_userEvent/blob/master/UserEvent20.cc) is the file to call in Phast. The number of the User Event is **20**. It calls the code of LCAnalysis.

## LCAnalysis

The main code. Performs the analysis on the mDSTs and outputs structured TTree with two principal branches: **DISEvents** and **Hadrons**. For each event (entry) of the tree corresponds a **DISEvent** and a list of **Hadrons**.

## MC treatment

The same framework can be used for MC. The outputs will contain reconstructed (similar to real data) and generated (MC) data.

## Usage

 - For Real Data :
```Bash
$DIR_PHAST/phast -u20 -c $OUTPUT_DIR -o $OUTPUT_DIR/phast_out.00.root -h $OUTPUT_DIR/hist_out.00.root -T data -U 00 -l [LIST OF FILES]
```

 - For Monte Carlo Data :
 ```Bash
$DIR_PHAST/phast -u20 -o $OUTPUT_DIR/phast_out.5125.root -h $OUTPUT_DIR/hist_out.5125.root -T MC2016 -U 5125 -l [LIST OF FILES]
 ```

## Further informations

 - [PHAST](http://ges.home.cern.ch/ges/phast/index.html)
