import math

# Step 1: Define constants for the problem
MATERIAL_W_CM = 12
MATERIAL_H_CM = 11

D1_DIAMETER_CM = 2
D1_CAPACITY_GB = 1

D2_DIAMETER_CM = 4
D2_CAPACITY_GB = 5

# Step 2: Calculate the optimal number of disks to maximize capacity
# We prioritize D2 disks as they offer more capacity per cm^2.
# We align the cuts with the 12cm and 11cm sides.
print("--- Calculating Maximum Storage Capacity ---")
print(f"Material dimensions: {MATERIAL_W_CM}cm x {MATERIAL_H_CM}cm")

# Calculate how many D2 disks (4x4cm) can fit
num_d2_along_w = MATERIAL_W_CM // D2_DIAMETER_CM
num_d2_along_h = MATERIAL_H_CM // D2_DIAMETER_CM
num_d2 = num_d2_along_w * num_d2_along_h
print(f"Number of D2 (4cm diameter) disks that fit: {num_d2_along_w} x {num_d2_along_h} = {num_d2} disks")

# The D2 disks occupy a 12cm x 8cm area
area_d2_w = num_d2_along_w * D2_DIAMETER_CM
area_d2_h = num_d2_along_h * D2_DIAMETER_CM
print(f"Area used by D2 disks: {area_d2_w}cm x {area_d2_h}cm")

# This leaves a remaining strip of material
remaining_strip_w = MATERIAL_W_CM
remaining_strip_h = MATERIAL_H_CM - area_d2_h
print(f"Remaining material is a {remaining_strip_w}cm x {remaining_strip_h}cm strip.")

# Now, calculate how many D1 disks (2x2cm) can fit in the remaining strip
num_d1_along_rem_w = remaining_strip_w // D1_DIAMETER_CM
num_d1_along_rem_h = remaining_strip_h // D1_DIAMETER_CM
num_d1 = num_d1_along_rem_w * num_d1_along_rem_h
print(f"Number of D1 (2cm diameter) disks from remaining material: {num_d1_along_rem_w} x {num_d1_along_rem_h} = {num_d1} disks")

# Calculate total capacity in Gigabytes (GB)
total_capacity_gb = (num_d2 * D2_CAPACITY_GB) + (num_d1 * D1_CAPACITY_GB)
print(f"\nTotal storage capacity = ({num_d2} D2 disks * {D2_CAPACITY_GB} GB) + ({num_d1} D1 disks * {D1_CAPACITY_GB} GB) = {total_capacity_gb} GB")

# Step 3: Calculate the size of a single data point
print("\n--- Calculating Data Point Size ---")
# Example record: "16 pm: extreme cold\n"
time_len = 5 # max length from "16 pm"
separator_len = 2 # ": "
category_len = 12 # max length from "extreme cold"
newline_len = 1 # for the newline character
data_point_size_bytes = time_len + separator_len + category_len + newline_len
print(f"A single data point record (e.g., '16 pm: extreme cold\\n') takes up a maximum of {data_point_size_bytes} bytes.")

# Step 4: Calculate the total number of data points that can be stored
print("\n--- Calculating Total Data Points ---")
BYTES_PER_GB = 10**9 # Using the standard that 1 GB = 1 billion bytes for storage
total_capacity_bytes = total_capacity_gb * BYTES_PER_GB
print(f"Total capacity in bytes = {total_capacity_gb} GB * {BYTES_PER_GB:,} bytes/GB = {total_capacity_bytes:,} bytes")

# Final calculation
max_data_points = total_capacity_bytes // data_point_size_bytes

print("\nThe final equation is:")
equation = f"( ( {num_d2} D2 disks * {D2_CAPACITY_GB} GB/disk ) + ( {num_d1} D1 disks * {D1_CAPACITY_GB} GB/disk ) ) * {BYTES_PER_GB:,} bytes/GB / {data_point_size_bytes} bytes/data point = {max_data_points:,} data points"
print(equation)

# Final answer formatted as requested
print("\nThe highest number of data points that can be collected is:")
print(f"{max_data_points:,}")

<<<1800000000>>>