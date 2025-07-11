import math

# Step 1: Define the initial parameters
material_width_cm = 12
material_height_cm = 11

d1_radius_cm = 1
d1_capacity_gb = 1
d2_radius_cm = 2
d2_capacity_gb = 5

d1_diameter_cm = d1_radius_cm * 2
d2_diameter_cm = d2_radius_cm * 2

# An estimation for the size of a single data point entry (e.g., "16 pm: extreme hot\n")
bytes_per_data_point = 25
# Standard conversion for storage capacity
gb_to_bytes = 10**9

print("--- Step 1: Maximizing storage by prioritizing D2 disks ---")
# Step 2: Calculate the number of D2 disks that can be cut from the material
# We use integer division to find how many disks fit along each dimension
num_d2_along_width = material_width_cm // d2_diameter_cm
num_d2_along_height = material_height_cm // d2_diameter_cm
total_d2_disks = num_d2_along_width * num_d2_along_height

print(f"Material dimensions: {material_width_cm}cm x {material_height_cm}cm")
print(f"D2 disk diameter: {d2_diameter_cm}cm")
print(f"Number of D2 disks that fit along the width ({material_width_cm}cm): {material_width_cm} // {d2_diameter_cm} = {num_d2_along_width}")
print(f"Number of D2 disks that fit along the height ({material_height_cm}cm): {material_height_cm} // {d2_diameter_cm} = {num_d2_along_height}")
print(f"Total D2 disks produced = {num_d2_along_width} * {num_d2_along_height} = {total_d2_disks}")
print("-" * 20)

# Calculate the capacity from D2 disks
capacity_from_d2_gb = total_d2_disks * d2_capacity_gb

# Step 3: Calculate remaining material and the number of D1 disks
# The D2 disks use a 12cm x 8cm area (3*4cm by 2*4cm)
area_used_by_d2_height = num_d2_along_height * d2_diameter_cm
remaining_material_width = material_width_cm
remaining_material_height = material_height_cm - area_used_by_d2_height

print("--- Step 2: Calculating D1 disks from remaining material ---")
print(f"The {total_d2_disks} D2 disks use an area of {material_width_cm}cm x {area_used_by_d2_height}cm.")
print(f"Remaining material is a strip of: {remaining_material_width}cm x {remaining_material_height}cm")

# Calculate how many D1 disks fit in the remaining strip
num_d1_along_width = remaining_material_width // d1_diameter_cm
num_d1_along_height = remaining_material_height // d1_diameter_cm
total_d1_disks = num_d1_along_width * num_d1_along_height

print(f"D1 disk diameter: {d1_diameter_cm}cm")
print(f"Number of D1 disks that fit along the width ({remaining_material_width}cm): {remaining_material_width} // {d1_diameter_cm} = {num_d1_along_width}")
print(f"Number of D1 disks that fit along the height ({remaining_material_height}cm): {remaining_material_height} // {d1_diameter_cm} = {num_d1_along_height}")
print(f"Total D1 disks produced = {num_d1_along_width} * {num_d1_along_height} = {total_d1_disks}")
print("-" * 20)

# Calculate the capacity from D1 disks
capacity_from_d1_gb = total_d1_disks * d1_capacity_gb

# Step 4: Calculate total capacity
print("--- Step 3: Calculating Total Storage Capacity ---")
total_capacity_gb = capacity_from_d2_gb + capacity_from_d1_gb
print(f"Capacity from D2 disks: {total_d2_disks} disks * {d2_capacity_gb} GB/disk = {capacity_from_d2_gb} GB")
print(f"Capacity from D1 disks: {total_d1_disks} disks * {d1_capacity_gb} GB/disk = {capacity_from_d1_gb} GB")
print(f"Total capacity = {capacity_from_d2_gb} GB + {capacity_from_d1_gb} GB = {total_capacity_gb} GB")
print("-" * 20)

# Convert total capacity to bytes
total_capacity_bytes = total_capacity_gb * gb_to_bytes

# Step 5: Calculate the total number of data points
print("--- Step 4: Calculating the Maximum Number of Data Points ---")
print(f"Total capacity in bytes = {total_capacity_gb} GB * {gb_to_bytes} bytes/GB = {total_capacity_bytes} bytes")
print(f"Estimated size of one data point: {bytes_per_data_point} bytes")

total_data_points = total_capacity_bytes // bytes_per_data_point

print("\nFinal Calculation:")
print(f"Highest number of data points = Total Capacity (bytes) / Bytes per Data Point")
print(f"Highest number of data points = {total_capacity_bytes} / {bytes_per_data_point}")
print(f"Result = {total_data_points}")
