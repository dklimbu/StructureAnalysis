# StructureAnalysis
Program to Calculate Structural Properties of Amorphous Materials
* *******************************************************************************  
 *   Program to calculate structural properties of amorphous silicon. The program reads a configuration of amorphous silicon model as a single XYZ file and calculates pair-correlation function (PCF), bond-angle distriution (BAD) and coordination number (CN).
 
 *   DIL LIMBU, USM
 *   APRIL 2018
 
 *   COMPILE:: 
 *             icc -o structure structure.cpp
 
 *   USAGE:: 
 *             ./structure INPUT_XYZ BIN_WIDTH
 
 *   OUTPUT :: 
 *             gr.dat  <-  pair correlation fuction;
 *             bad.dat <-  bond-angle districution;
 *             CN.dat  <-  coordination number 
 * *******************************************************************************  

<p>
  <img src="gr.png" width="400" height="300" align=left>
  <img src="bad.png" width="400" height="300" align=left>
</p>
