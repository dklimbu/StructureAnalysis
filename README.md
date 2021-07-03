# StructureAnalysis
Program to Calculate Structural Properties of Amorphous Materials. 
*****************************************************************
This program calculates structural properties of amorphous silicon
by reading a configuration of amorphous silicon model as 
XYZ file and calculates Pair-Correlation Function (PCF)
Structure Factor (S(k)), Bond-Angle Distribution (BAD), and 
generate Coordination Number (CN) and Nearest Neighor Map (NMAP)

* DIL LIMBU, USM
* APRIL 2018

* TO COMPILE:: Use Makefile provided in /src dir

* USAGE:: ./structure.x INPUT.XYZ

* OUTPUTS:: 
*           gr.dat   <- Pair Correlation Function
*           Sk.dat   <- Structure Factor
*           bad.dat  <- Bond Angle Distribution
*           Cn.dat   <- Coordination Number
*           nmap.dat <- List of ALL Nearest Neighbors
*************************************************************
* PLOT::
*      RUN Python script to generate plots
*************************************************************
<p>
  <img src="gr.png" width="400" height="300" align=left>
   <img src="Sk.png" width="400" height="300" align=left>
  <img src="bad.png" width="400" height="300" align=left>
</p>
