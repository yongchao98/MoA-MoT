import math

# Define material and disk properties
material_w_cm = 12
material_h_cm = 11

d1_diameter_cm = 2
d1_capacity_gb = 1

d2_diameter_cm = 4
d2_capacity_gb = 5

# Estimate the size of a single data point record in bytes
# e.g., "12 pm: extreme hot\n" is less than 20 bytes.
data_point_size_bytes = 20

# Step 1: Maximize the number of D2 disks (the more capacity-dense type)
# Place D2 disks (4cm diameter) on the 12cm x 11cm material
num_d2_along_w = material_w_cm // d2_diameter_cm
num_d2_along_h = material_h_cm // d2_diameter_cm
total_d2_disks = num_d2_along_w * num_d2_along_h

# Step 2: Calculate remaining material and fit D1 disks
# The D2 disks use an area of (3 * 4cm) x (2 * 4cm) = 12cm x 8cm
used_h_cm = num_d2_along_h * d2_diameter_cm
rem_material_w_cm = material_w_cm
rem_material_h_cm = material_h_cm - used_h_cm

# Now, fit D1 disks (2cm diameter) into the remaining 12cm x 3cm strip
num_d1_along_rem_w = rem_material_w_cm // d1_diameter_cm
num_d1_along_rem_h = rem_material_h_cm // d1_diameter_cm
total_d1_disks = num_d1_along_rem_w * num_d1_along_rem_h

# Step 3: Calculate the total storage capacity in GB
total_capacity_gb = (total_d2_disks * d2_capacity_gb) + (total_d1_disks * d1_capacity_gb)

# Step 4: Convert total capacity to bytes
bytes_per_gb = 1024**3
total_capacity_bytes = total_capacity_gb * bytes_per_gb

# Step 5: Calculate the total number of data points that can be recorded
max_data_points = total_capacity_bytes // data_point_size_bytes

# --- Output the results ---
print(f"Optimal packing strategy:")
print(f"- Number of D2 disks (5GB): {total_d2_disks}")
print(f"- Number of D1 disks (1GB): {total_d1_disks}")
print(f"This gives a total storage capacity of {total_capacity_gb} GB.")
print("\nCalculating the maximum number of data points:")
print(f"({total_capacity_gb} GB * {bytes_per_gb} bytes/GB) / {data_point_size_bytes} bytes/data_point = {max_data_points} data points")

print(f"\nThe highest number of data points that can be collected and recorded is {max_data_points}.")
print("Final equation with all numbers:")
print(f"{total_capacity_gb} * {bytes_per_gb} / {data_point_size_bytes} = {max_data_points}")