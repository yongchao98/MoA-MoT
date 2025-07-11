import math

# Step 1: Define the initial parameters for the material and disks.
material_width_cm = 12
material_height_cm = 11

d1_diameter_cm = 2
d1_capacity_gb = 1

d2_diameter_cm = 4
d2_capacity_gb = 5

# Step 2: Calculate the optimal number of disks to maximize capacity.
# We prioritize the larger D2 disks first.
num_d2_along_width = material_width_cm // d2_diameter_cm
num_d2_along_height = material_height_cm // d2_diameter_cm
total_d2_disks = num_d2_along_width * num_d2_along_height
capacity_from_d2_gb = total_d2_disks * d2_capacity_gb

# Calculate the area remaining after placing the D2 disks.
width_used_by_d2 = num_d2_along_width * d2_diameter_cm
height_used_by_d2 = num_d2_along_height * d2_diameter_cm

remaining_height_cm = material_height_cm - height_used_by_d2
# The remaining area is a strip of material_width_cm x remaining_height_cm

# Now, fit the smaller D1 disks into the remaining strip.
num_d1_in_remaining_width = material_width_cm // d1_diameter_cm
num_d1_in_remaining_height = remaining_height_cm // d1_diameter_cm
total_d1_disks = num_d1_in_remaining_width * num_d1_in_remaining_height
capacity_from_d1_gb = total_d1_disks * d1_capacity_gb

# Calculate the total maximum capacity.
total_capacity_gb = capacity_from_d2_gb + capacity_from_d1_gb

# Step 3: Define data point size and conversion factor.
data_point_size_bytes = 32
bytes_per_gb = 1024 * 1024 * 1024

# Step 4: Calculate the total capacity in bytes.
total_capacity_bytes = total_capacity_gb * bytes_per_gb

# Step 5: Calculate the final number of data points.
max_data_points = total_capacity_bytes // data_point_size_bytes

# Step 6: Print the detailed, step-by-step solution.
print("--- Part 1: Calculating Maximum Storage Capacity ---")
print(f"Calculating the number of D2 disks (4cm diameter) from a 12cm x 11cm sheet:")
print(f"Fit along width: {material_width_cm} // {d2_diameter_cm} = {num_d2_along_width} disks")
print(f"Fit along height: {material_height_cm} // {d2_diameter_cm} = {num_d2_along_height} disks")
print(f"Total D2 disks = {num_d2_along_width} * {num_d2_along_height} = {total_d2_disks}")
print(f"Capacity from D2 disks = {total_d2_disks} * {d2_capacity_gb} GB = {capacity_from_d2_gb} GB")
print("-" * 20)
print(f"The D2 disks use a {width_used_by_d2}cm x {height_used_by_d2}cm area, leaving a {material_width_cm}cm x {remaining_height_cm}cm strip.")
print(f"Calculating the number of D1 disks (2cm diameter) from the remaining strip:")
print(f"Fit along width: {material_width_cm} // {d1_diameter_cm} = {num_d1_in_remaining_width} disks")
print(f"Fit along height: {remaining_height_cm} // {d1_diameter_cm} = {num_d1_in_remaining_height} disk")
print(f"Total D1 disks = {num_d1_in_remaining_width} * {num_d1_in_remaining_height} = {total_d1_disks}")
print(f"Capacity from D1 disks = {total_d1_disks} * {d1_capacity_gb} GB = {capacity_from_d1_gb} GB")
print("-" * 20)
print(f"Total Maximum Capacity = {capacity_from_d2_gb} GB (from D2) + {capacity_from_d1_gb} GB (from D1) = {total_capacity_gb} GB")

print("\n--- Part 2: Calculating Total Number of Data Points ---")
print(f"Total capacity in bytes = {total_capacity_gb} GB * {bytes_per_gb} bytes/GB = {total_capacity_bytes} bytes")
print(f"Estimated size of a single data point = {data_point_size_bytes} bytes")
print("\nFinal Equation:")
print(f"Max Data Points = Total Capacity (bytes) / Size per Data Point (bytes)")
print(f"                = {total_capacity_bytes} / {data_point_size_bytes}")
print(f"                = {max_data_points}")