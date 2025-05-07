from paraview.simple import *
import glob
import os
import csv
import sys
import yaml

yml_file = sys.argv[1]
vtu_file = sys.argv[2]
output_file = sys.argv[3]

# ========== Get parameters from .yml file ==========
with open(yml_file, 'r') as f:
    config = yaml.safe_load(f)

# Example: extract values from the config
h = config["grid"]["aligned_nano_slit_with_unstructured_reservoirs_2d"]["reservoir_height"]
r = config["grid"]["aligned_nano_slit_with_unstructured_reservoirs_2d"]["pore_radius"]
l = config["grid"]["aligned_nano_slit_with_unstructured_reservoirs_2d"]["pore_length"]
L_REF = config["non_dim"]["L_REF"]
h_non_dim = h/L_REF
r_non_dim = r/L_REF
l_non_dim = l/L_REF

# Disable automatic rendering to improve performance
paraview.simple._DisableFirstRenderCameraReset()

# Parameters to customize
vtu_file_path = vtu_file  # Path to your .pvd file
line_point1 = [l_non_dim*1.5, h_non_dim-r_non_dim, 0]              # Starting point of line
line_point2 = [l_non_dim*1.5, h_non_dim, 0]                # Ending point of line
num_samples = 1000                        # Number of points along the line

print(f"Reading PVD file: {vtu_file_path}")

# Check if file exists
if not os.path.exists(vtu_file_path):
    print(f"Error: File {vtu_file_path} not found!")
    exit(1)

# === Use XMLPartitionedUnstructuredGridReader for .pvt ===
reader = XMLUnstructuredGridReader(FileName=[vtu_file_path])
reader.UpdatePipeline()


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
csv_file = os.path.join("test.csv")
print(array_names)
# Write CSV without separate Time column
with open(csv_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    header = ["X", "Y", "Z"] + array_names
    writer.writerow(header)
    writer.writerows(rows)

print("Processing completed successfully")