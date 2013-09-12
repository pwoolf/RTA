Regular Task Analysis (RTA) 
===

This software is designed for identifying classification patterns in time series data.  Patterns are defined based on a limited set of regular expression skeletons.  These skeletons are exaustively searched and scored based on a Fisher's exact test.  The final scores are corrected using the Bonferoni correction to identify possible overfitting.  

Quickstart
===

Downalod, and run the synthetic case:

     python RegularTaskAnalysis.py -f synthetic_data.txt

Note that the code is parallelized, so will take advantage of all of the cores presen on the machine.  Note also that the code may take quite a while to run, depending on your dataset and machine.   



