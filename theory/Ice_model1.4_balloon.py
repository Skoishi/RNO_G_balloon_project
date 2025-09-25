from NuRadioMC.utilities import medium
from NuRadioMC.SignalProp.radioproparaytracing import radiopropa_ray_tracing
from AntPosCal.time_differences.travel_time_simulator_radiopropa import TravelTimeSimulatorRadioPropa
import mattak.Dataset  # RNO-G data handling
import time
import numpy as np
import argparse
import uproot  # ROOT file reading
import matplotlib.pyplot as plt
import gpxpy  # GPS track parsing
import gpxpy.gpx
import gzip
import pandas as pd
import os
import readRNOGDataMattak 
import logging
import requests
from NuRadioReco.detector import detector
from rnog_analysis_tools.coordinate_system.coordinate_system import CoordinateSystem
import datetime

t0 = time.time()

det = detector.Detector(source="rnog_mongo")
det.update(datetime.datetime(2024, 6, 1))
#station_deployed=[11,12,13,21,22,23,24]
station=21
station_dis=det.get_absolute_position(station)
ch_arr=[]
coordsys=CoordinateSystem()

# Instead of using pulser location from detector,
# use user-defined coordinate (Lat, Lon, Height)
user_x = -417.927690   # meters
user_y = 458.540055   # meters
user_height = 337.32  # meters

# Convert geodetic to ENU (relative to station position)
station_geodetic = coordsys.enu_to_geodetic(station_dis[0], station_dis[1], station_dis[2])

balloon_dis = [user_x, user_y, user_height]

print(f"Using user-defined coordinate: ENU = {balloon_dis}")


for ch in range(24):
    ch_dis=station_dis + det.get_relative_position(station,ch)
    print(ch_dis)
    ch_arr.append(np.array(ch_dis))
    #coord=coordsys.enu_to_geodetic(ch_dis[0],ch_dis[1],ch_dis[2])
    #ch_coord.append(np.array(coord))



det1 = detector.Detector(json_filename="/users/PAS2608/youwei/RNO_G/NuRadioMC/NuRadioReco/detector/RNO_G/RNO_season_2024.json")
det1.update(datetime.datetime(2024, 9, 14))

print(balloon_dis)

ch_index=np.array(range(24))


#-------------------------------------------------------------------------------------------------------------------
ice_simp = medium.greenland_simple()
ice_poly = medium.greenland_poly5()
ice_simp_propa=TravelTimeSimulatorRadioPropa(ice_simp)
ice_poly_propa=TravelTimeSimulatorRadioPropa(ice_poly)

print("data successfully loaded in")
#-------------------------------------------------------------------------------------------------------------------
# Using NumPy zeros
matrix_poly = pd.DataFrame(np.zeros((24, 24)))

sphere_sizes_i = np.array([250,25, 2., .5])

for i in ch_index:
    for j in ch_index:
        if i > j:
            try:

                
                # Calculate travel time for channel i

                #ice_poly_propa.set_sphere_sizes(sphere_sizes_i)
                time_i = ice_poly_propa.get_travel_time(ch_arr[i], balloon_dis)
                
                # Calculate travel time for channel j
                #ice_poly_propa.set_sphere_sizes(sphere_sizes_j)
                time_j = ice_poly_propa.get_travel_time(ch_arr[j], balloon_dis)
                
                # Compute time difference
                time_diff = time_i - time_j
                matrix_poly.iloc[i, j] = time_diff
                
                # Optional: Print for debugging
                print(f"Channels ({i},{j}): Time difference = {time_diff}")
                #print(f"  Used sphere sizes: {sphere_sizes_i} (ch{i}), {sphere_sizes_j} (ch{j})")
                
            except Exception as e:
                print(f"Error for ({i},{j}): {e}")
                matrix_poly.iloc[i, j] = np.nan


matrix_simp = pd.DataFrame(np.zeros((24, 24)))
print("matrix_simp successfully calculated")


for i in ch_index:
    for j in ch_index:
        if i > j:
            try:
                
                # Calculate travel time for channel i
                # ice_simp_propa.set_sphere_sizes(sphere_sizes_i)
                time_i = ice_simp_propa.get_travel_time(ch_arr[i], balloon_dis)
                
                # Calculate travel time for channel j
                #ice_simp_propa.set_sphere_sizes(sphere_sizes_j)
                time_j = ice_simp_propa.get_travel_time(ch_arr[j], balloon_dis)
                
                # Compute time difference
                time_diff = time_i - time_j
                matrix_simp.iloc[i, j] = time_diff
                
                # Optional: Print for debugging
                #print(f"Channels ({i},{j}): Time difference = {time_diff}")
                #print(f"  Used sphere sizes: {sphere_sizes_i} (ch{i}), {sphere_sizes_j} (ch{j})")
                
            except Exception as e:
                print(f"Error for ({i},{j}): {e}")
                matrix_poly.iloc[i, j] = np.nan
print("matrix_poly successfully calculated")

                



#-------------------------------------------------------------------------------------------------------------------
# Write matrices to file (with station number in filename)
with open(f'Ice_model/matrix_ice_model_poly_station_{station}_run_2205.txt', 'w') as f:
    f.write("Time delay matrix:\n")
    f.write(matrix_poly.to_string())

with open(f'Ice_model/matrix_ice_model_simp_station_{station}_run_2205.txt', 'w') as f:
    f.write("Time delay matrix:\n")
    f.write(matrix_simp.to_string())


t_total = time.time() - t0
print(t_total)