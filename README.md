# seq_processing
This repository contains a collection of Python-based tools developed by Chi-Chun Ho, initially intended for internal use at Division of Chemical Pathology, Department of Clinical Pathology, Pamela Youde Nethersole Eastern Hospital (PYNEH), Hong Kong, China.
The scripts are provided “as is”, without being actively supported or maintained.
1.  Installation
  * The scripts are written in Python 3 and depend on a number of modules. The easiest way to install the required modules is by installing Anaconda. Anaconda may be downloaded from https://www.anaconda.com/distribution/
  * Alternatively, the required modules may be individually installed (see below).
  * During installation, the Windows user is suggested to select installation of “Just Me” (as recommended), “Add Anaconda to my PATH environment variable” (despite not recommended) and “Register Anaconda as my default Python 3.x”.
  * Most scripts require a working Internet connection to function properly. For example, the core program autoprimer.py downloads transcript information from the Ensembl REST server, and other programs (e.g. web_search.py) may require other web services. 
2. Running
  * Check that Python is properly installed. In the Windows command processor (cmd.exe), type
```
python
```
and the Anaconda Python interpreter should be invoked:
```
Python 3.7.0 (default, Aug 14 2018, 19:12:50) [MSC v.1900 32 bit (Intel)] :: Anaconda, Inc. on win32
Type "help", "copyright", "credits" or "license" for more information.
>>>
```
You may then try to run the individual scripts. Running the scripts with no arguments provided displays a (hopefully informative) help message:
```
C:\Users\COMPUTER_USER\Documents\Scripts\seq_processing>python build_bed_file.py
usage: build_bed_file.py [-h] [-f CDS_FLANK] transcripts [transcripts ...]
build_bed_file.py: error: the following arguments are required: transcripts
# AutoPrimer program by Tom C.C. Ho (c) 2017
# Version 1.0j build 20190110
# Release note: Updated the Gene() class to allow automatic transcript selection, and bug fixes
```
Note that, in Linux, ```python``` should probably be replace by ```python3``` or similar.

3.   Additional dependencies
  * Some tools require additional Python modules
  * For batch_tag_exons.py and exon_coverage_report.py
```
python -m pip install biopython
python -m pip install pdfkit
python -m pip install numpy
python -m pip install matplotlib
```