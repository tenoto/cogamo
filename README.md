# cogamo package of the Thundercloud Project 

[![hackmd-github-sync-badge](https://hackmd.io/-2miC95hTvC5lLqzmqyZxw/badge)](https://hackmd.io/-2miC95hTvC5lLqzmqyZxw)

Analysis codes and examples for the Compact Gamma-ray Monitor (Coamo) of the "Thundercloud Project" organized by the GROWTH Collaboration. 

## 1. What is the Thundercloud Project? 

The Thundercloud Project is a scientific program to observe high-energy phenomena occurring in thunderclouds using Cogamo, a gamma-ray measuring instrument installed on the ground. The program has been performed in Kanazawa, Japan, to observe thunderclouds that arrive every winter. This project employs a form of Citizen Science, and here we share the basic analysis source codes with supporters. For details of the project, see the paper (Enoto et al., PTEP, 2022, in prep). A story of the project was featured by Nature, ["Mystery gamma rays could help solve age-old lightning puzzle"](https://www.nature.com/articles/d41586-021-00395-3). The publications of the GROWTH collaboration are compiled in [the ADS Library](https://ui.adsabs.harvard.edu/public-libraries/NEXUZ5uATYaFPJQxcvs9lQ). 

## 2. Installation and setup

For your first installation, use the git clone command.
```
git clone https://github.com/tenoto/cogamo.git
```

Then, edit the directory paths in the "setenv/setenv.bashrc" file, and read this file under the top directory of the "cogamo" library,
```
source setenv/setenv.bashrc
```

## 3. Raw data format

All the measured data are recorded in a MicroSD card of the Cogamo detector. 

### 3-1. Directory structure of the Cogamo micro SD
The cogamo data are organized in the order of detector IDs (Det_ID). In the following case, the Det_ID 11 has "config.csv", "data", and "log" directories. 
```
011
├── config.csv
├── data
│   ├── 011_20200305_00.csv
│   ├── 011_20200305_01.csv
...
│   └── 011_20200305_23.csv
└── log
    ├── 011_20181002.csv
    ├── 011_20181031.csv
...
```
The "data" folder includes the event data, while the "log" folder include the house keeping (HK) data.

### 3-2. Event data 
The event data file recorded all the individual radiation events. The file name is 

```
[det_id]_[yyyymmdd]_[hour].csv
```

where the time is defined in JST.

The event data file include the following columns:



1. minute
2. sec
3. 1/10000 sec
4. ADC channel (pha: pulse height amplitude) [0-1023]


### 3-3. House Keeping (HK) data 
The House Keeping (HK) data file includes basic parameters of a detector and its environment. The file name is 

```
[det_id]_[yyyymmdd].csv
```

where the time is defined in JST. The file is generated per a day.

The raw csv-format file is included the following columns (',' separated):

1. yyyy-mm-dd (JST)
2. HH:MM:SS (JST)
3. data recording interval (min)
4. rate1 (cps) below "AREABD1" of " config.csv"
5. rate2 (cps) between "AREABD1" and  "AREABD2" of " config.csv"
6. rate3 (cps) between "AREABD2" and  "AREABD3" of " config.csv"
7. rate4 (cps) between "AREABD3" and  "AREABD4" of " config.csv"
8. rate5 (cps) between "AREABD4" and  "AREABD5" of " config.csv"
9. rate6 (cps) above "AREABD5" of " config.csv"
10. temperature (degC)
11. pressure (hPa)
12. humidity (%)
13. the maximum value among the difference of 10-sec moving average of count rates to the latest count rates (10秒移動平均とCPS値との差の最大値:定義を確認する必要がある) 
14. optical illumination (lux)
15. gps status (0:invalid, 1 or 2: valid)
16. longitude (deg)
17. latitude (deg)

### 3-4. Config file
The config file is incldued in a MicroSD card of the Cogamo detector, which setup the instrumental parameters

- INTERVAL,5 (*Internal for writing of the HK data*)
- ID,100 (*Detector ID*)
- AREABD1,1600 (*Threshold energy between rate1 and rate2 in a keV unit (i.e., 1600 keV)*)
- AREABD2,3000 (*Threshold energy between rate2 and rate3 in a keV unit (i.e., 3000 keV)*)
- AREABD3,5000 
- AREABD4,8000 
- AREABD5,10000 
- TIMECONS,10 (*Not used, same in the following parameters*)
- SPCTRINT,120 
- MODE,0 
- NBIN,60 
- MULTIP,1.5 
- NPRE,200 
- RTH,300 
- RBACK,280

## 4. Functions of the main library

To be prepared

## 5. Examples 

To be prepared

## 6. Cogamo response file

In X-rays and gamma-ray bands, individual photons deposit their total energy to the detector (total absorption) and sometimes leave part of the incident energy (e.g., Compton scattering). Thus, without considering such a complicated detector response, we cannot determine the incident energy spectrum and flux before interacting with the detector. This problem is solved by preparing a two-dimensional matrix, called the "detector response," in astrophysics and working to find the model that best reproduces the measured gamma-ray spectrum by fitting.

Here, we perform data analysis using the spectral analysis software Xspec, used in X-ray and gamma-ray astronomy, a part of the HEASoft package developed by NASA.

## 6-1. Response generation [You can skip this section]

Most of you can skip this subsection if you do not intend to generate your own response. This library already includes prepared response files generated as written below. 

The Thundercloud Project uses Geant4 simulations to generate a gamma-ray photon list of the detector interaction. This list has at least three columns, 'event number,' 'incident energy (keV),' and 'deposited energy (keV) to the scintillator.' 

In the setup file (see setenv.bashrc, this gamma-ray photon list is placed in the $COGAMO_RESPONSE_DATA_PATH directory. For example, 'CsI_ext_source_coated.txt' is calculated for the Cogamo FY2021 model (5x5x15 cm3 CsI coated by the Al package). This includes 10^7 photons with random initial energy in 0.04-41 MeV. 

The Xspec-format response file is generated in the following command ('cogamo/respnose' directory).
```
generate_cogamo_response.py parameter/cogamo_fy2020_flat_rsp.yaml
```
where 'parameter/cogamo_fy2020_flat_rsp.yaml' describes input parameters. This 'generate_cogamo_response.py' script fills the gamma-ray photon list into the two-dimensional matrix in the response file and further incorporates the energy degradation as a function of the photon energy. This command makes 'cogamo_fy2020_flat.rsp' (flat injection into the cogamo). This file is also stored in the 'cogamo/response' directory. 

You can try this in the test field. 

```
cd test/response    
script/generate_cogamo_respons_flat.py
```

References (for you to learn more about Xspec response files)
- [Creating FITS format files for XSPEC](http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/fits/fitsfiles.html)
- [The Calibration Requirements for Spectral Analysis](http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3.2)

## 6-2. Spectral fitting example 

Example spectral files are generated via 'test/single_evtdata/script/run.py' command from the 016_20210108_06 data. These example files are placed under 'cogamo/example/xspec' directory. 


```
%> xspec
XSPEC12>data 1 016_20210108_06_bst01_src_bin.pha 
XSPEC12>back 1 016_20210108_06_bst01_bgd.pha 
XSPEC12>resp 1 ../../cogamo/response/cogamo_fy2020_flat.rsp 
XSPEC12>setplot energy mev
XSPEC12>ignore **-0.2,20.0-**
XSPEC12>cpd /xw
XSPEC12>pl ld
XSPEC12>model cutoffpl
1:cutoffpl:PhoIndex>1.4
2:cutoffpl:HighECut>8000,0.1,0.1,0.1,10000,10000
3:cutoffpl:norm>10
XSPEC12>renorm
XSPEC12>query yes
XSPEC12>fit
XSPEC12>pl ld del
XSPEC12>flux 200 20000
```
This gives flux of 3.0e-6 ergs/cm2/s in the 0.2-20 MeV. 


