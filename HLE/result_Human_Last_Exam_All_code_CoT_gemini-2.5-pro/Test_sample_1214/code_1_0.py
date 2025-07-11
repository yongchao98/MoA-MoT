import math

# Step 1: Define constants for the problem
material_w = 12  # cm
material_h = 11  # cm

d1_radius = 1  # cm
d1_capacity_gb = 1  # GB
d1_diameter = d1_radius * 2

d2_radius = 2  # cm
d2_capacity_gb = 5  # GB
d2_diameter = d2_radius * 2

# Step 2: Estimate the size of one data point
# Example: "16 pm: extreme hot" + newline character is approx. 20 bytes
data_point_size_bytes = 20
bytes_per_gb = 10**9

# Step 3: Find the optimal number of disks to maximize capacity
# We place the larger D2 disks first in a grid pattern
num_d2_along_w = math.floor(material_w / d2_diameter)
num_d2_along_h = math.floor(material_h / d2_diameter)
num_d2_disks = num_d2_along_w * num_d2_along_h

# The D2 disks occupy a rectangular area. The rest of the material is a strip.
# D2s take up (3 * 4cm) by (2 * 4cm) = 12cm by 8cm.
# The remaining material is a strip of 12cm by (11cm - 8cm) = 12cm by 3cm.
rem_material_w = material_w
rem_material_h = material_h - (num_d2_along_h * d2_diameter)

# Now, calculate how many D1 disks can fit in this remaining strip
num_d1_along_w = math.floor(rem_material_w / d1_diameter)
num_d1_along_h = math.floor(rem_material_h / d1_diameter)
num_d1_disks = num_d1_along_w * num_d1_along_h

# Step 4: Calculate the total storage capacity
total_capacity_gb = (num_d2_disks * d2_capacity_gb) + (num_d1_disks * d1_capacity_gb)
total_capacity_bytes = total_capacity_gb * bytes_per_gb

# Step 5: Calculate the maximum number of data points
max_data_points = total_capacity_bytes / data_point_size_bytes

# Step 6: Print the results clearly, showing the final equation
print("To maximize storage, the best combination is:")
print(f"- {num_d2_disks} D2 disks (5GB each)")
print(f"- {num_d1_disks} D1 disks (1GB each) from the remaining material")
print("-" * 30)
print("This gives the following total capacity:")
print(f"Total GB = ({num_d2_disks} * {d2_capacity_gb}) + ({num_d1_disks} * {d1_capacity_gb}) = {total_capacity_gb} GB")
print("-" * 30)
print("Calculating the maximum number of data points:")
print(f"Total Bytes = {total_capacity_gb} GB * {bytes_per_gb} Bytes/GB = {total_capacity_bytes} Bytes")
print(f"Size per Data Point = {data_point_size_bytes} Bytes")
print("-" * 30)
print("Final Equation:")
print(f"Maximum Data Points = ({total_capacity_bytes}) / ({data_point_size_bytes}) = {int(max_data_points)}")
