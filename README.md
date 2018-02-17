# AUC calculation
<b>This script calculates the area-under-the-curve of user-defined gene sets for a given dataset </b>  
<b>Author</b>: Mudra Hegde  
<b>Email</b>: mhegde@broadinstitute.org  
<b>Version: 1.1 </b>  
  
<b>Required packages</b>
1. pandas
2. numpy
3. scipy
4. statsmodels
5. scikit-learn
  
<b>Inputs</b>
1. <b>Input File</b>: .txt file with list of sgRNAs in the first column and log-fold changes for every sample in the following columns. Column names of log-fold changes will be sample name in the output file. 
2. <b>Chip File</b>: .txt file to map sgRNAs to gene symbols;First column should be the list of sgRNAs and second column should be the corresponding gene symbols
3. <b>Gene sets folder</b>: Folder with .txt files of each gene set for AUC calculation. Column name for the gene set in the output file will be the file name of the gene set. Gene set file should be a .txt file with one column of gene symbols and no header.
4. <b>Output File</b>
  
<b>To run this script, type the following on the terminal:</b>  
python calc_auc_v1.1.py --input-file \<Path to inputfile\> --chip-file \<Path to chip file\> --gene-set-folder \<Path to folder\> --outputfile \<Path to output file\>  
    
To view the available options in the script, type the following on the terminal:  
python calc_auc_v1.1.py -h 

Sets of essential and non-essential genes from Hart et al.,2015 and Hart et al.,2014 are included along with this script.  


  

