import math

# This script calculates the maximum number of data points that can be stored on HDDs
# made from a 12x11cm rectangular sheet of material.

# --- Step 1: Define problem constants ---

# Material dimensions in cm
MATERIAL_W = 12
MATERIAL_H = 11

# Disk properties
# D1: 1cm radius (2cm diameter), 1GB capacity
D1_DIAMETER = 2
D1_CAPACITY_GB = 1
# D2: 2cm radius (4cm diameter), 5GB capacity
D2_DIAMETER = 4
D2_CAPACITY_GB = 5

# Assumption: 1 GB = 1,000,000,000 bytes (common for HDD manufacturers)
BYTES_PER_GB = 1_000_000_000

# --- Step 2: Estimate the size of a single data point ---

# A data point is a string like "12 am: extreme cold" plus a newline character.
# The longest possible category is "extreme cold" (12 chars).
# The time format can be "12 am" (5 chars).
# The full string is "12 am: " (7 chars) + "extreme cold" (12 chars) + "\n" (1 char).
# Total size = 7 + 12 + 1 = 20. Let's use the provided example for more accuracy.
# "12 am: extreme cold\n" has a length of 22 characters. We'll use this as our size in bytes.
DATA_POINT_SIZE_BYTES = len("12 am: extreme cold\n")


# --- Step 3: Find the optimal disk combination for maximum capacity ---

def find_max_capacity(rect_w, rect_h):
    """
    Calculates max capacity by trying all row combinations of D2 disks
    and filling the rest of the space with D1 disks.
    """
    max_capacity = 0
    best_config = {'d2': 0, 'd1': 0}

    # Iterate through the number of possible rows of D2 disks
    for num_d2_rows in range(math.floor(rect_h / D2_DIAMETER) + 1):
        # Calculate how many D2 disks fit
        d2s_per_row = math.floor(rect_w / D2_DIAMETER)
        num_d2s = num_d2_rows * d2s_per_row
        
        # Calculate the remaining rectangular area for D1 disks
        remaining_h = rect_h - (num_d2_rows * D2_DIAMETER)
        
        # Calculate how many D1 disks fit in the remaining area
        d1s_per_row = math.floor(rect_w / D1_DIAMETER)
        num_d1_rows = math.floor(remaining_h / D1_DIAMETER)
        num_d1s = d1s_per_row * num_d1_rows
        
        # Calculate the total capacity for this configuration
        current_capacity = (num_d2s * D2_CAPACITY_GB) + (num_d1s * D1_CAPACITY_GB)
        
        # If this configuration is better, save it
        if current_capacity > max_capacity:
            max_capacity = current_capacity
            best_config = {'d2': num_d2s, 'd1': num_d1s}
            
    return max_capacity, best_config

# To be thorough, check both orientations of the material (12x11 and 11x12)
cap1, config1 = find_max_capacity(MATERIAL_W, MATERIAL_H)
cap2, config2 = find_max_capacity(MATERIAL_H, MATERIAL_W)

# Select the best result from the two orientations
if cap1 > cap2:
    max_total_capacity_gb = cap1
    best_num_d2 = config1['d2']
    best_num_d1 = config1['d1']
else:
    max_total_capacity_gb = cap2
    best_num_d2 = config2['d2']
    best_num_d1 = config2['d1']


# --- Step 4: Calculate the total number of data points and print the result ---

# Total capacity in bytes
total_storage_bytes = max_total_capacity_gb * BYTES_PER_GB
# Total number of data points (using integer division)
total_data_points = total_storage_bytes // DATA_POINT_SIZE_BYTES

# Print the final calculation as a single, comprehensive equation
print("The final calculation to find the highest number of data points is:")
print(f"((Number of D2 disks * Capacity of D2) + (Number of D1 disks * Capacity of D1)) * Bytes per GB / Bytes per data point")
print(f"(({best_num_d2} * {D2_CAPACITY_GB}) + ({best_num_d1} * {D1_CAPACITY_GB})) * {BYTES_PER_GB:,} / {DATA_POINT_SIZE_BYTES} = {total_data_points:,}")
