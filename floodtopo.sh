#!/usr/bin/env bash

module load CDO netCDF-Fortran

basetopo=/work/kaplan_lab/datasets/topography/GEBCO2025/GEBCO_globe.nc

icefile=/work/kaplan_lab/datasets/ICE-7G/I7G_NA.VM7_1deg.6.nc 

# interpolate ICE-7G topographic anomaly to target grid 

cdo -f nc4 -P 8 remapbic,$basetopo -selname,stgit,Topo_Diff $icefile tmp.nc
