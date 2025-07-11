import math

# Step 1: Define the constants for the problem.
# Material dimensions in cm
material_width = 12
material_height = 11

# Disk D1 (small) properties
d1_diameter = 2  # cm, requires a 2x2 cm square
d1_capacity_gb = 1  # GB

# Disk D2 (large) properties
d2_diameter = 4  # cm, requires a 4x4 cm square
d2_capacity_gb = 5  # GB

# Data storage constants
# Assumed size for an efficient data point (8-byte timestamp + 1-byte category ID)
data_point_size_bytes = 9
# Conversion factor for GB to bytes (1024^3)
bytes_per_gb = 1073741824

# Step 2: Calculate the optimal number of each disk type to maximize capacity.
# We prioritize D2 disks because they have a higher capacity density (5GB/16cm^2 > 1GB/4cm^2).
# We first calculate how many D2 disks (4x4cm squares) can fit.
num_d2_along_width = math.floor(material_width / d2_diameter)
num_d2_along_height = math.floor(material_height / d2_diameter)
num_d2 = num_d2_along_width * num_d2_along_height

# Calculate the area used by D2 disks and the remaining area.
height_used_by_d2 = num_d2_along_height * d2_diameter
remaining_height = material_height - height_used_by_d2
remaining_width = material_width

# Step 3: Calculate how many D1 disks (2x2cm squares) fit in the remaining area.
num_d1_along_width = math.floor(remaining_width / d1_diameter)
num_d1_along_height = math.floor(remaining_height / d1_diameter)
num_d1 = num_d1_along_width * num_d1_along_height

# Step 4: Calculate the total storage capacity in GB and then in bytes.
total_capacity_gb = (num_d2 * d2_capacity_gb) + (num_d1 * d1_capacity_gb)
total_capacity_bytes = total_capacity_gb * bytes_per_gb

# Step 5: Calculate the maximum number of data points.
max_data_points = total_capacity_bytes // data_point_size_bytes

# Step 6: Print the calculation step-by-step, showing the numbers in the final equation.
print(f"To maximize capacity, we calculated the optimal mix of disks.")
print(f"Number of D2 disks (5GB): {num_d2}")
print(f"Number of D1 disks (1GB) in remaining area: {num_d1}")
print("-" * 30)

print("Final Calculation:")
# First, show the equation for total Gigabytes
print(f"Total GB = ({num_d2} * {d2_capacity_gb}) + ({num_d1} * {d1_capacity_gb}) = {total_capacity_gb} GB")
# Next, show the equation for total data points
print(f"Max Data Points = ({total_capacity_gb} GB * {bytes_per_gb} bytes/GB) / {data_point_size_bytes} bytes/point")
# Finally, show the full equation with the result
print(f"({total_capacity_gb} * {bytes_per_gb}) // {data_point_size_bytes} = {max_data_points}")
