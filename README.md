# cogamo
Analysis package for the Compact Gamma-ray Monitor (Coamo) of the "Thundercloud Project" organized by the GROWTH Collaboration.

This page summarizes example of the Cogamo data analyses. 

## Installation and setup

For your first installation, use the git clone command.
```
git clone https://github.com/tenoto/cogamo.git
```

Then, under the top directory of  the "cogamo" library,
```
source setenv/setenv.bashrc
```

## Raw data format

All the measured data are recorded in a MicroSD card of the Cogamo detector. 

### Directory structure
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

### House Keeping (HK) data 
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


### Event data 
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

### Config file
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

