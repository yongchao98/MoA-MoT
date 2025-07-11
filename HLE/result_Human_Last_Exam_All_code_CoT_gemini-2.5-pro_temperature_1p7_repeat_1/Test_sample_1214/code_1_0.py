import math

# Step 1: Determine the optimal disk configuration and maximum capacity.
# Based on geometric packing analysis, using D2 disks is optimal.
# A 12x11cm sheet can fit 8 D2 disks (4cm diameter).

num_disks = 8
capacity_per_disk_gb = 5 # GB

# Total capacity in Gigabytes (GB)
total_capacity_gb = num_disks * capacity_per_disk_gb

# Step 2: Define the size of a single data point in bytes.
# We assume a fixed size for each record to account for variations and overhead.
data_point_size_bytes = 32 # bytes

# Step 3: Calculate the total number of data points.
# Convert total capacity to bytes. 1 GB = 1024 * 1024 * 1024 bytes.
bytes_in_a_gb = 1024**3
total_capacity_bytes = total_capacity_gb * bytes_in_a_gb

# Calculate the maximum number of data points
max_data_points = total_capacity_bytes // data_point_size_bytes

# Print the final result showing the full equation
print("This plan maximizes storage by creating 8 D2 disks from the material.")
print(f"Number of disks = {num_disks}")
print(f"Capacity per disk = {capacity_per_disk_gb} GB")
print(f"Assumed size per data point = {data_point_size_bytes} bytes")
print("-" * 20)
print("The highest number of data points can be calculated as:")
print(f"({num_disks} disks * {capacity_per_disk_gb} GB/disk * {bytes_in_a_gb} bytes/GB) / {data_point_size_bytes} bytes/point")
print(f"= {max_data_points}")
