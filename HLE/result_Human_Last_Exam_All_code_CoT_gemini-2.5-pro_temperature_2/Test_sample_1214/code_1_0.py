import math

# Step 1: Define problem constants
material_w_cm = 12
material_h_cm = 11

d1_radius_cm = 1
d1_capacity_gb = 1

d2_radius_cm = 2
d2_capacity_gb = 5

# A single data point like "16 pm: extreme hot" is a string. We'll assume a fixed size
# for each record, which is a reasonable and safe estimate for this kind of data.
data_point_size_bytes = 32

# Step 2: Determine the square footprint for each disk
d1_side = d1_radius_cm * 2
d2_side = d2_radius_cm * 2

# Step 3: Calculate the number of disks based on the packing strategy
# Prioritize D2 disks as they are more capacity-dense.
# We consider both orientations of the material sheet to ensure we find the maximum.

# Function to calculate layout for a given orientation
def calculate_layout(width, height):
    # Fit D2 disks first
    num_d2_w = math.floor(width / d2_side)
    num_d2_h = math.floor(height / d2_side)
    num_d2 = num_d2_w * num_d2_h
    
    # Calculate remaining area. The used area is (num_d2_w * d2_side) x (num_d2_h * d2_side)
    # The remaining area is composed of two rectangles.
    # Rectangle 1:
    rem_w1 = width
    rem_h1 = height - (num_d2_h * d2_side)
    num_d1_in_rem1 = math.floor(rem_w1 / d1_side) * math.floor(rem_h1 / d1_side)
    
    # Rectangle 2:
    rem_w2 = width - (num_d2_w * d2_side)
    rem_h2 = num_d2_h * d2_side
    num_d1_in_rem2 = math.floor(rem_w2 / d1_side) * math.floor(rem_h2 / d1_side)

    total_d1 = num_d1_in_rem1 + num_d1_in_rem2
    
    total_capacity = (num_d2 * d2_capacity_gb) + (total_d1 * d1_capacity_gb)
    
    return num_d2, total_d1, total_capacity

# Calculate for both orientations
num_d2_orient1, num_d1_orient1, capacity1 = calculate_layout(material_w_cm, material_h_cm)
num_d2_orient2, num_d1_orient2, capacity2 = calculate_layout(material_h_cm, material_w_cm)

# Choose the orientation with the highest capacity
if capacity1 > capacity2:
    num_d2 = num_d2_orient1
    num_d1 = num_d1_orient1
    total_capacity_gb = capacity1
else:
    num_d2 = num_d2_orient2
    num_d1 = num_d1_orient2
    total_capacity_gb = capacity2


# Step 4: Calculate total storage capacity in bytes
bytes_per_gb = 1024**3
total_capacity_bytes = total_capacity_gb * bytes_per_gb

# Step 5: Calculate maximum number of data points
max_data_points = total_capacity_bytes // data_point_size_bytes

# Print the results step-by-step
print(f"--- HDD Manufacturing Calculation ---")
print(f"Material size: {material_w_cm}cm x {material_h_cm}cm")
print(f"D1 disk (radius {d1_radius_cm}cm) requires a {d1_side}x{d1_side}cm square and stores {d1_capacity_gb}GB.")
print(f"D2 disk (radius {d2_radius_cm}cm) requires a {d2_side}x{d2_side}cm square and stores {d2_capacity_gb}GB.")
print("\n--- Optimal Packing Result ---")
print(f"Number of D2 disks produced: {num_d2}")
print(f"Number of D1 disks produced: {num_d1}")

print("\n--- Storage Capacity Calculation ---")
print(f"Total GB from D2 disks: {num_d2} * {d2_capacity_gb} = {num_d2 * d2_capacity_gb} GB")
print(f"Total GB from D1 disks: {num_d1} * {d1_capacity_gb} = {num_d1 * d1_capacity_gb} GB")
print(f"Total storage capacity: {total_capacity_gb} GB")
print(f"Total storage in bytes: {total_capacity_gb} GB * {bytes_per_gb} Bytes/GB = {total_capacity_bytes} Bytes")

print("\n--- Data Point Calculation ---")
print(f"Assumed size of a single data point record: {data_point_size_bytes} Bytes")
print("Final equation:")
print(f"(({num_d2} * {d2_capacity_gb}) + ({num_d1} * {d1_capacity_gb})) * {bytes_per_gb} / {data_point_size_bytes} = {max_data_points}")
print(f"\nThe highest number of data points that can be collected is: {max_data_points}")
