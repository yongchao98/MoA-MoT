import math

# Step 1: Define material and disk properties
rect_l = 12  # cm
rect_w = 11  # cm

d1_diameter = 2  # cm
d1_capacity = 1  # GB

d2_diameter = 4  # cm
d2_capacity = 5  # GB

# Step 2: Calculate the optimal number of disks using the mixed packing strategy
# First, pack the more space-efficient D2 disks
num_d2_disks = math.floor(rect_l / d2_diameter) * math.floor(rect_w / d2_diameter)

# Calculate remaining area dimensions
rem_l = rect_l
rem_w = rect_w - (math.floor(rect_w / d2_diameter) * d2_diameter)

# Second, pack D1 disks in the remaining rectangular area
num_d1_disks = math.floor(rem_l / d1_diameter) * math.floor(rem_w / d1_diameter)

# Step 3: Calculate the maximum total storage capacity
max_capacity_gb = (num_d2_disks * d2_capacity) + (num_d1_disks * d1_capacity)

# Step 4: Define the size of a single data point record
# Based on the longest example "12 am: extreme cold" plus a newline character
bytes_per_datapoint = 20  # bytes

# Step 5: Calculate the final number of data points
# Using the standard conversion: 1 GB = 1024^3 bytes
bytes_per_gb = 1024**3
total_datapoints = (max_capacity_gb * bytes_per_gb) // bytes_per_datapoint

# --- Output of the final calculation ---
print("The highest number of data points can be calculated as follows:")
print(f"\n1. Maximize storage capacity by mixing disk types:")
print(f"   - Number of D2 disks: {num_d2_disks}")
print(f"   - Number of D1 disks in remaining space: {num_d1_disks}")
print(f"   - Maximum Capacity = ({num_d2_disks} * {d2_capacity}) + ({num_d1_disks} * {d1_capacity}) = {max_capacity_gb} GB")

print(f"\n2. Calculate total data points:")
print(f"   - Total capacity in bytes = {max_capacity_gb} GB * {bytes_per_gb} bytes/GB")
print(f"   - Size of one data point record = {bytes_per_datapoint} bytes")

print(f"\nFinal Equation:")
print(f"Total Points = ({max_capacity_gb} * {bytes_per_gb}) / {bytes_per_datapoint} = {total_datapoints}")

<<<1932735283>>>