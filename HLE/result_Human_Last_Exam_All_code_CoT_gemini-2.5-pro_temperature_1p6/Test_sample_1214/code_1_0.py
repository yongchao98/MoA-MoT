import math

# Step 1: Maximize Storage Capacity
# Material and Disk Dimensions
material_w = 12  # cm
material_h = 11  # cm

d1_radius = 1  # cm
d1_diameter = 2 * d1_radius
d1_capacity = 1  # GB

d2_radius = 2  # cm
d2_diameter = 2 * d2_radius
d2_capacity = 5  # GB

# To maximize capacity, we should prioritize the disk with the higher capacity-per-area.
# D1 occupies a 2x2=4 cm^2 square, gives 1 GB. Value = 1/4 = 0.25 GB/cm^2.
# D2 occupies a 4x4=16 cm^2 square, gives 5 GB. Value = 5/16 = 0.3125 GB/cm^2.
# D2 is more space-efficient, so we place as many D2 disks as possible first.

# Calculate how many D2 disks can fit in a grid layout.
# Orientation 1:
d2_along_w = material_w // d2_diameter
d2_along_h = material_h // d2_diameter
num_d2_config1 = d2_along_w * d2_along_h

# The best grid arrangement for D2 disks (4cm diameter) on a 12x11 cm sheet is
# a 3x2 grid, occupying a 12cm x 8cm area.
num_d2 = 3 * 2
print("Step 1: Calculate the maximum storage capacity from the material.")
print(f"Based on the 12cm x 11cm material, we can create an optimal mix of disks.")
print(f"Maximizing the more space-efficient D2 disks, we can fit {num_d2} D2 disks.")

# This arrangement (3 D2 disks along the 12cm side, 2 along the 11cm side)
# occupies a 12cm x 8cm area.
width_used_by_d2 = d2_along_w * d2_diameter
height_used_by_d2 = d2_along_h * d2_diameter
print(f"This uses a {width_used_by_d2}cm x {height_used_by_d2}cm area of the material.")

# The remaining area is a 12cm x 3cm strip.
rem_w = material_w
rem_h = material_h - height_used_by_d2
print(f"The remaining area is a {rem_w}cm x {rem_h}cm strip.")

# Now, we fill the remaining area with D1 disks (2cm diameter).
d1_in_rem_w = rem_w // d1_diameter
d1_in_rem_h = rem_h // d1_diameter
num_d1 = d1_in_rem_w * d1_in_rem_h
print(f"In this remaining area, we can fit {num_d1} D1 disks.")

# Calculate the total capacity
total_capacity_gb = (num_d2 * d2_capacity) + (num_d1 * d1_capacity)
print("\nCalculating total capacity:")
print(f"Total Capacity = ({num_d2} D2 disks * {d2_capacity} GB/disk) + ({num_d1} D1 disks * {d1_capacity} GB/disk)")
print(f"Total Capacity = {num_d2 * d2_capacity} GB + {num_d1 * d1_capacity} GB = {total_capacity_gb} GB")


# Step 2: Calculate Data Point Size
# A data point consists of: Time (e.g., "12 pm", 5 chars) + Separator (": ", 2 chars) + Category ("extreme cold", 12 chars) + Newline (1 char)
size_time = 5
size_separator = 2
size_category = 12  # Longest category is "extreme cold"
size_newline = 1
data_point_size_bytes = size_time + size_separator + size_category + size_newline

print("\nStep 2: Calculate the size of a single data point.")
print(f"The longest data point string is '12 pm: extreme cold\\n', which is {data_point_size_bytes} characters long.")
print(f"Assuming 1 character = 1 byte, one data point requires {data_point_size_bytes} bytes.")

# Step 3: Calculate Total Number of Data Points
# Standard conversion for HDD capacity is 1 GB = 1,000,000,000 bytes.
gb_to_bytes = 10**9
total_capacity_bytes = total_capacity_gb * gb_to_bytes
num_data_points = total_capacity_bytes // data_point_size_bytes

print("\nStep 3: Calculate the total number of data points that can be stored.")
print("The final calculation is (Total Capacity in Bytes) / (Size per Data Point).")
# Printing the full equation as requested
print(f"Equation: (({num_d2} * {d2_capacity}) + ({num_d1} * {d1_capacity})) * {gb_to_bytes} / {data_point_size_bytes}")
print(f"Result = {num_data_points:,}") # Using a comma for readability in this printout
print("\nHighest number of data points that can be collected and recorded:")
print(num_data_points)
