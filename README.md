# cogamo
Analysis package of CoGaMo detectors for the Thundercloud Project 


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
