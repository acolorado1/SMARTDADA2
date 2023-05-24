# SMARTDADA2

- [SMARTDADA2](#smartdada2)
  - [Repo Directory](#repo-directory)
  - [Installations and Dependencies](#installations-and-dependencies)
  - [Workflow](#workflow)
    - [To Run Program](#to-run-program)
    - [Input](#input)
    - [Output](#output)
  - [Contact](#contact)

## Repo Directory

|Directory Name|Description|
|---------------|-----------|
|.github| Contains github actions tests|
|notebooks| Contains series of steps and explanation of our method|
|smartdada2| Contains main scripts that will allow user to better understand the quality of their reads|

## Installations and Dependencies

To install this program:

1. clone repository
2. create environment
3. activate environment
4. install smartdada2 module into python environment

For example in the terminal type:

```python
git clone https://github.com/acolorado1/DADA2ParameterExploration.git
conda env create -f smartdada2_env.yaml
conda activate {environemtnname i.e., smartdada2}
pip install -e .
```

## Workflow

### To Run Program

This program can be run using snakemake. In the Snakefile you must put the file path of the FASTQ file of your choosing and adjust parameters as needed. Parameters include:  

- PAIRED (bool): Are you using paired or signle end reads?
- SUBSAMPLE (int): Number of reads you want to use when creating visualizations. 
- AVG_Q_SCORE (default = 30.0): Determines when obvious trimming will stop.
- OBV_TRIMMING_MAX (default = 0.1): Determines when obvious trimming will throw a warning that the read might be too short. Default warns if trimming is over 10% of the read on either end.
- MAX_TRIMMING (default = 0.2): Calculates max index trimming on either end when finding average expected error for position. Default will only calculate up to 20% trim/truncating on each end.

<img width="682" alt="Screen Shot 2023-05-23 at 1 58 10 PM" src="https://github.com/acolorado1/SMARTDADA2/assets/68305443/6ca259fb-35d6-4498-9828-65984619aff8">

Once that is done write in the terminal:

```
snakemake -c 1 
```

*Note:* ```-c``` specifies the max number of cores you want to use.


### Input

This script takes in one FASTQ formatted file. For example:

```
@M07186:25:000000000-DHMKP:1:1101:15331:1350 1:N:0:1
TNCGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCGCGTCGCTTGTGAAAGCCCGGGGCTTAACTCCGGGTCTGCAGGCGATACGGGCATAACTTGAGTGCTGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGGAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCAGTAACTGACGCTGAGGAGCGAAAGCATGGGGAGCGAAC
+
>#>>>>>AA1>AEEEGEEGGFDFEGGEGHGHHHHGGGGGGGHHHFHHGGGGGFFFEFHGHBEGGCEE>EGGHHBBF1GGHGGGGGGHGHFHDHGGCBCGCGFH0EC/CGCGCGCGHEFCGHHGC=<GGDDDGHGG?-.EGHACCHFGFFFBCF.BB9?@;-AEBFFFGG?@-FFF9/;/BEBFFAFFF@@@FF--@@=?EFF<FFFFFFEF?-;FFFBF//FBBAB-FEBAE==@;B9BFFFFBA-->@@-
@M07186:25:000000000-DHMKP:1:1101:15094:1354 1:N:0:1
TNCGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACA
+
A#>>>ABBBAABFGGGFGGGGGFDGGGGHFHHGHHHFFGGHGHFEEGECGGGGGGFGGEGFHHHHHGG@EGFGHFHHHHHGGGGGGHHFGHHHHGHHGHFHHHHHHHEE@GHHHHGGGHGHFGHHHHGHFHHHHGGGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFACDFAFFDAFFFFFFFFFFFFFFFFFFDDDDFFFFBFDDCFFFFFFFFFFFFFFFF
@M07186:25:000000000-DHMKP:1:1101:14910:1354 1:N:0:1
TNCGATTAACCCAAACTAATAGGCCTACGGCGTAAAGCGTGTTCAAGATACTTTTACACTAAAGTTAAAACTTAACTAAGCCGTAAAAAGCTACAGTTATCATAAAATAAACCACGAAAGTGACTTTATAATAATCTGACTACACGATAGCTAAGACCCAAACTGGGATTAGAAACCCCTGTAGTCCGGCTGGCTGACTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAATAGACGTGCTAGGTAT
+
A#>>>AABFFFB2AAEGGGGGGGHCGFHFG?EAEEAFGGGFEGH5AGGHFHHHHHHFFHHHGGHHHHHHGHHHHBDGGFHHHGGGGGGEDHHHHGHHGHHHHHFGEFGGFFFFFEFGGEF?FDEFHHHGHFFHGGHHEFFHHGHHDGEGGHGFHHGHHHHHHHBFGCGHHHHFDGGHFGFGCDGHHHDCGGCCABEGGHHH0CGHGGG/FGGG@FFGGGGFEFFB0CFBB<-;--.////.....//.9//
```

### Output

An interactive output containing resulting images and two TSV files are output from these scripts in a directory called *forward/output* and *reverse/output*. In addition, static formats of the images will also be present in the subdirectory called **plots**. 

The interactive output will appear as an html file that will resemble the following: 

https://github.com/acolorado1/SMARTDADA2/assets/68305443/9f44536e-c4d5-466d-85b6-fca4a6c92a02

The two TSVs that will be output and that will be used to create the figures look like the following:

TrimInfo:

```
LeftIndex	RightIndex	ReadLength	AvgEEPerPosition	      RightTrim	LeftTrim
0	        201	        201	        0.00025149061996092216	234	      0
0	        202	        202	        0.0002525037876296975	  234	      0
0	        203	        203	        0.0002531623874528494	  234	      0
0	        204	        204	        0.00025375036528248743	234	      0
```

SumEEInfo:

```
NoTrimming         ObviousTrimming
1.7478821923032126 1.545527496734793
0.686050226996564   0.6824597001220604
1.6827705639478587 0.9463769578124168
```

## Contact

Angela Sofia Burkhart Colorado - angelasofia.burkhartcolorado@cuanschutz.edu

Erik Serrano - erik.serrano@cuanschutz.edu
