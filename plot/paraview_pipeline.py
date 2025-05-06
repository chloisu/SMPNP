from paraview.simple import *
import glob
import os
import csv

# Disable automatic rendering to improve performance
paraview.simple._DisableFirstRenderCameraReset()

# Parameters to customize
pvd_file_path = '../../../../Data/20250417_r_dependency_h80/r1/out.pvd'  # Path to your .pvd file
line_point1 = [0, 0, 0]              # Starting point of line
line_point2 = [1, 0, 0]                # Ending point of line
num_samples = 100                          # Number of points along the line
output_dir = '.'   # Path to save the CSV file

print(f"Reading PVD file: {pvd_file_path}")

# Check if file exists
if not os.path.exists(pvd_file_path):
    print(f"Error: File {pvd_file_path} not found!")
    exit(1)

# Set up the reader for PVD files
reader = PVDReader(registrationName='PVDReader', FileName=pvd_file_path)

# Update the pipeline to load data
reader.UpdatePipeline()

# Get available timesteps
try:
    timesteps = reader.TimestepValues
    if not timesteps:
        print("No timesteps found in the PVD file")
        exit(1)
    
    # Get the last timestep
    last_timestep = timesteps[-1]
    print(f"Found {len(timesteps)} timesteps. Using last timestep: {last_timestep}")
    
except AttributeError:
    print("No timesteps found, processing as single timestep")
    last_timestep = 0.0

# Update reader to the last timestep
reader.UpdatePipeline(last_timestep)

# Sample data along the line
plot_line = PlotOverLine(
    Input=reader,
    Point1=line_point1,
    Point2=line_point2,
    Resolution=num_samples
)
plot_line.UpdatePipeline()

# Fetch VTK table
table = servermanager.Fetch(plot_line)

array_names = []  # will be populated after first extraction
# Determine quantity names once
if not array_names:
    pt_data = table.GetPointData()
    array_names = [pt_data.GetArray(i).GetName() for i in range(pt_data.GetNumberOfArrays())]

# Collect rows for this timestep (omit Time column since it is in filename)
rows = []
for i in range(table.GetNumberOfPoints()):
    x, y, z = table.GetPoints().GetPoint(i)
    # Only spatial coords and field values
    row = [x, y, z]
    for name in array_names:
        arr = table.GetPointData().GetArray(name)
        row.append(arr.GetValue(i))
    rows.append(row)

# Define output CSV per timestep, embedding time in filename
# Example: profile_t0.1234.csv
csv_file = os.path.join(output_dir, "test.csv")

# Write CSV without separate Time column
with open(csv_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    header = ["X", "Y", "Z"] + array_names
    writer.writerow(header)
    writer.writerows(rows)

print("Processing completed successfully")