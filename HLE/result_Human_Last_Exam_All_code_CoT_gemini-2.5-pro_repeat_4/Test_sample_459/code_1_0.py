import math

# --- Given Parameters ---
num_levels = 6
largest_level_size_gb = 1
mem_buffer_size_kb = 1
insert_rate_bps = 16000  # bytes per second
page_size_bytes = 2500

# --- Convert units to bytes for consistency ---
largest_level_size_bytes = largest_level_size_gb * (1024**3)
mem_buffer_size_bytes = mem_buffer_size_kb * 1024

print("--- Parameters ---")
print(f"Number of levels (L): {num_levels}")
print(f"Largest level size: {largest_level_size_gb} GB")
print(f"Memory buffer size: {mem_buffer_size_kb} KB")
print(f"Insert rate: {insert_rate_bps} bytes/s")
print(f"Page size: {page_size_bytes} bytes")
print("-" * 20 + "\n")


# --- Step 1: Calculate the size ratio (T) ---
# In a leveled LSM tree, Size(L_i) = T * Size(L_{i-1}).
# Assuming 6 levels on disk (L0 to L5), we have Size(L5) = T^5 * Size(L0).
# We assume the effective size of the base level (L0) is based on the memory buffer size.
print("--- Step 1: Calculate Size Ratio (T) ---")
# The number of multiplicative steps from the base level to the largest level is (num_levels - 1)
power = num_levels - 1
size_ratio_t = (largest_level_size_bytes / mem_buffer_size_bytes) ** (1/power)
# Round to the nearest integer as T is typically an integer value.
size_ratio_t = round(size_ratio_t)

print(f"The size relationship is: Largest_Level_Size = T^(Num_Levels - 1) * Base_Level_Size")
print(f"T = ({largest_level_size_bytes} / {mem_buffer_size_bytes}) ^ (1 / {power})")
print(f"Calculated size ratio (T): {size_ratio_t}")
print("-" * 20 + "\n")


# --- Step 2: Calculate Total I/O Rate in Bytes/s ---
# For leveled compaction, the I/O for moving data from L_{i-1} to L_i is:
# Read from L_{i-1}: insert_rate
# Read from L_i: T * insert_rate
# Write to L_i: (T+1) * insert_rate
# Total I/O for one merge: (1 + T + T + 1) * insert_rate = (2T + 2) * insert_rate
# This happens for (num_levels - 1) merge operations.
# The total I/O is the initial write to L0 plus the I/O from all merges.
print("--- Step 2: Calculate Total I/O Rate (Bytes/s) ---")
num_merges = num_levels - 1
io_per_merge_factor = (2 * size_ratio_t + 2)
total_io_bytes_per_sec = insert_rate_bps + num_merges * io_per_merge_factor * insert_rate_bps

print("I/O for one merge = (Read_prev_lvl + Read_curr_lvl + Write_curr_lvl)")
print(f"I/O for one merge = (1 + {size_ratio_t} + ({size_ratio_t}+1)) * Insert_Rate = {io_per_merge_factor} * Insert_Rate")
print(f"Total I/O = Initial_Write + Num_Merges * I/O_per_Merge")
print(f"Total I/O Rate = {insert_rate_bps} + {num_merges} * {io_per_merge_factor} * {insert_rate_bps} = {total_io_bytes_per_sec} bytes/s")
print("-" * 20 + "\n")


# --- Step 3: Calculate Total Page I/O Rate ---
# Convert the total byte I/O rate to page I/O rate.
print("--- Step 3: Calculate Minimum Total Page I/O Rate ---")
page_io_rate = total_io_bytes_per_sec / page_size_bytes
print("Page I/O Rate = Total I/O Rate (bytes/s) / Page Size (bytes)")
print(f"The equation with all values is:")
print(f"({insert_rate_bps} + ({num_levels} - 1) * (2 * {size_ratio_t} + 2) * {insert_rate_bps}) / {page_size_bytes}")
print(f"\nMinimum total page I/O rate: {page_io_rate:.1f} pages/s")
