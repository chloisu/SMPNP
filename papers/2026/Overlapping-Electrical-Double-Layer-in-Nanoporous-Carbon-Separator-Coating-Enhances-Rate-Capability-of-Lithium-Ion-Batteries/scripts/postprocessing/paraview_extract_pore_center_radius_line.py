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

start_x = l_non_dim*1.5
end_x = l_non_dim*1.5
start_y = h_non_dim-r_non_dim
end_y = h_non_dim
# Disable automatic rendering to improve performance
paraview.simple._DisableFirstRenderCameraReset()

# Parameters to customize
vtu_file_path = vtu_file  # Path to your .pvd file
line_point1 = [start_x, start_y+1e-6, 0]              # Starting point of line
line_point2 = [end_x, end_y-1e-6, 0]                # Ending point of line
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

pt_data = table.GetPointData()
array_names = []
array_components = []

for i in range(pt_data.GetNumberOfArrays()):
    arr = pt_data.GetArray(i)
    name = arr.GetName()
    n_comp = arr.GetNumberOfComponents()
    if n_comp == 1:
        array_names.append(name)
        array_components.append((name, None))  # Scalar field
    else:
        for c in range(n_comp):
            comp_name = f"{name}:{c}"
            array_names.append(comp_name)
            array_components.append((name, c))  # Vector component

# Collect rows for this timestep (omit Time column since it is in filename)
rows = []
for i in range(table.GetNumberOfPoints()):
    x, y, z = table.GetPoints().GetPoint(i)
    row = [x, y, z]
    for base_name, comp in array_components:
        arr = pt_data.GetArray(base_name)
        if comp is None:
            value = arr.GetValue(i)
        else:
            value = arr.GetComponent(i, comp)
        row.append(value)
    rows.append(row)

# Define output CSV per timestep, embedding time in filename
# Example: profile_t0.1234.csv
# Write CSV without separate Time column
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    header = ["X", "Y", "Z"] + array_names
    writer.writerow(header)
    writer.writerows(rows)

print("Processing completed successfully")