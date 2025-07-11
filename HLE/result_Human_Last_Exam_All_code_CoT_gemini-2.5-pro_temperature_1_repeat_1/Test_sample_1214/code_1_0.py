import math

# --- Plan ---
# 1. Define material and disk specifications.
# 2. Determine the optimal number of D1 and D2 disks to maximize capacity from the material.
#    - This involves solving a 2D packing problem. We prioritize the more space-efficient D2 disks.
#    - We check both 12x11 and 11x12 orientations of the material to find the true maximum.
# 3. Calculate the total storage capacity in GB and then convert it to bytes.
# 4. Estimate the size of a single data point record.
# 5. Divide the total byte capacity by the single data point size to get the final answer.

# --- Step 1: Define Constants ---

# Material dimensions in cm
material_w = 12
material_h = 11

# Disk D1 properties
d1_diameter = 2  # cm (1cm radius * 2)
d1_capacity_gb = 1  # GB

# Disk D2 properties
d2_diameter = 4  # cm (2cm radius * 2)
d2_capacity_gb = 5  # GB

# Data point size estimation in bytes
# "16 pm: extreme hot\n" -> time (8 bytes) + category (16 bytes) + separator (1 byte)
data_point_size_bytes = 25

# Conversion factor for GB to Bytes
bytes_per_gb = 1024**3

# --- Step 2: Maximize Storage Capacity ---

# We use a function to calculate the best layout for a given material orientation
def calculate_best_layout(mat_w, mat_h):
    """
    Calculates the number of D1 and D2 disks by prioritizing the denser D2 disks first,
    then filling the remaining space with D1 disks.
    """
    # Pack as many D2 (4x4) disks as possible
    num_d2_along_w = mat_w // d2_diameter
    num_d2_along_h = mat_h // d2_diameter
    num_d2 = num_d2_along_w * num_d2_along_h

    # The D2 disks occupy a rectangular area
    d2_area_w = num_d2_along_w * d2_diameter
    d2_area_h = num_d2_along_h * d2_diameter

    # The remaining area is an L-shape, which we split into two rectangles to fill with D1 disks
    
    # Rectangle 1 (along the width of the D2 area)
    rem_h1 = mat_h - d2_area_h
    num_d1_rect1 = (d2_area_w // d1_diameter) * (rem_h1 // d1_diameter)
    
    # Rectangle 2 (the rest of the material's width)
    rem_w2 = mat_w - d2_area_w
    num_d1_rect2 = (rem_w2 // d1_diameter) * (mat_h // d1_diameter)

    num_d1 = num_d1_rect1 + num_d1_rect2
    
    return num_d2, num_d1

# The packing result is the same for 12x11 and 11x12, so we calculate it once.
num_d2_disks, num_d1_disks = calculate_best_layout(material_w, material_h)

# --- Step 3: Calculate Total Capacity ---

total_capacity_gb = (num_d1_disks * d1_capacity_gb) + (num_d2_disks * d2_capacity_gb)
total_capacity_bytes = total_capacity_gb * bytes_per_gb

# --- Step 4 & 5: Calculate and Print Final Result ---

# Calculate the maximum number of data points (cannot be a fraction)
max_data_points = total_capacity_bytes // data_point_size_bytes

print("--- Step-by-Step Calculation ---")
print(f"1. Maximizing Storage from a {material_w}x{material_h} cm Sheet:")
print(f"   - Number of D2 disks (5GB, 4cm diameter) that can be made: {num_d2_disks}")
print(f"   - Number of D1 disks (1GB, 2cm diameter) from remaining material: {num_d1_disks}")
print("")
print("2. Calculating Total Capacity:")
print(f"   - Total GB = ({num_d1_disks} D1 disks * {d1_capacity_gb} GB) + ({num_d2_disks} D2 disks * {d2_capacity_gb} GB)")
print(f"   - Final Equation for GB: ({num_d1_disks} * {d1_capacity_gb}) + ({num_d2_disks} * {d2_capacity_gb}) = {total_capacity_gb} GB")
print("")
print("3. Converting to Bytes:")
print(f"   - Total Bytes = {total_capacity_gb} GB * {bytes_per_gb} Bytes/GB")
print(f"   - Final Equation for Bytes: {total_capacity_gb} * {bytes_per_gb} = {total_capacity_bytes} Bytes")
print("")
print("4. Calculating Maximum Data Points:")
print(f"   - Estimated size of one data point record: {data_point_size_bytes} Bytes")
print(f"   - Total Data Points = Total Bytes / Bytes per Point")
print(f"   - Final Equation: {total_capacity_bytes} / {data_point_size_bytes} = {max_data_points} data points")
print("")
print(f"The highest number of data points that can be collected is: {max_data_points}")
