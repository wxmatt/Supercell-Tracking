# Supercell-Tracking
Scripts used to track supercell within the Rasmussen and Liu 2016 dataset.

Mainclass.py does the heavy lifting of handling the reflectivity, UH, Bunkers motion, etc. and then passing that information into the TiNT tracking.
This is performed on a daily basis (1200 UTC to 1200 UTC). The final product for is daily .csv files of all convective cell tracks within the domain for that day. 
All convective cells are tracked that fall within the TiNT parameters (so no explicit supercell classification is performed here).

PyartGridder.py simply converts the WRF reflectivity into a TiNT usable format. 

Combine.py connects and does some bookkeeping for each daily track .csv file produced by Mainclass.py. The final product from this script is a .csv file and dictionary
that contains all supercells. This script includes the supercell classification. 
