# SCALE-LETKF

**SCALE-LETKF** is a data assimilation and numerical weather prediction system utilizing [LETKF](https://github.com/takemasa-miyoshi/letkf). 

## Sources 

### SCALE model

[Scalable Computing for Advanced Library and Environment (SCALE)](https://scale.riken.jp/) - Regional Model (RM)

* [Get source code](https://scale.riken.jp/download/)  
  The current SCALE-LETKF version supports SCALE-5.4.5.

* [Documentation](http://scale.riken.jp/doc/)  
  Model description and user guide.

* [Data files (source of topo/landuse)](https://scale.riken.jp/download/#datasets)

### Database

Input data for testcase is located at the following places.

**Fugaku:** /share/hp150019/scale_database  
**Hibuna:** /home/amemiya/scale_database  

(**For fugaku users**: Be sure to copy the database to your directory under /vol0003 or /vol0004, as /share is not accesible directly from the compute node.)   

## [Getting started](Getting-started.md) 

## Run testcases

### Ideal experiments

* [Baroclinic waves](Baroclinic-waves.md)

* [Supercell](Supercell.md)

### Real experiments

* [PREPBUFR DA for Japan area](18km_Japan)

* [PAWR DA](PAWR-DA.md)

* [Himawari 8 radiances DA](Himawari-8-radiances-DA.md)

* [PAWR DA with the merged version](PAWR-DA-with-the-merged-version.md)

## Run your experiments

## [Configuable settings](Configuable-settings.md)

## Others

### [Realtime system](Realtime-system.md)

### [Old documents](Old-documents.md)
