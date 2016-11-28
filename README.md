# repeat_finding
# This program is to find all possible segments that may be associated with gene deletion and insertion.

Files and folders in this project are listed as follows:

0. this readme file;
1. Folder for blast software which is applied in this project to make alignment between two gene blocks;
2. Folder for genome strain data which includes a set of reference genome;
3. faNamesPre.txt with genome name information;
4. cordis_disp_abstract.txt. This file is an abstract for maf file. You'd better exchange it with a complete maf file for processing;
5. executable file sort_result to process the blast alignment result. 
6. extractSq.sh file to extract a segment from a reference with the information of start position and end position;
7. 1_step.cpp for the first step of finding repeat. The usage of this file is written at the beginning of this file;
8. 2_step.cpp for the second step of fiding repeat. The usage of this file is written at the beginning of this file.

After processing, a folder called ./compare_result will be generated. Delete all the empty folders and files inside it and get the final folder that store all the possible repeats that may be related to gene indels.
