import math

# --- Given Parameters ---
num_levels = 6
largest_level_size_gb = 1
mem_buffer_size_kb = 1
insert_rate_bps = 16000  # bytes per second
page_size_bytes = 2500

# --- Step 1: Determine LSM Tree Geometry ---

# Convert all sizes to bytes for consistent calculations
mem_buffer_bytes = mem_buffer_size_kb * 1024
largest_level_bytes = largest_level_size_gb * (1024**3)

# The number of levels on disk is num_levels - 1 (L1 to L5)
# The number of scaling steps from memory (L0) to the largest level (L5) is num_levels - 1
num_scaling_steps = num_levels - 1

# Calculate the size ratio 'T'
# largest_level_size = mem_buffer_size * (T ^ num_scaling_steps)
size_ratio_T = (largest_level_bytes / mem_buffer_bytes) ** (1 / num_scaling_steps)
# Since (1024^2)^(1/5) = (2^20)^(1/5) = 2^4 = 16, T should be an integer.
# We'll round it to be safe, but it calculates to 16.0
size_ratio_T = round(size_ratio_T)

print(f"Calculated size ratio T between levels: {size_ratio_T}")
print("-" * 30)

# --- Step 2: Model the I/O Rate ---

# I/O for the L0 (memory) to L1 (disk) flush is purely write I/O
# The average rate of this I/O is the insert rate.
io_rate_mem_flush_bps = insert_rate_bps

# The number of merges between disk levels is (num_levels - 2)
# These are L1->L2, L2->L3, L3->L4, L4->L5
num_disk_to_disk_merges = num_levels - 2

# For each disk-to-disk merge, the average I/O rate is 2 * (T+1) * R_insert
# This accounts for reading both levels and writing the new merged level.
# I/O Rate = (read_L_i + read_L_i+1 + write_L_i+1) / (time_to_fill_L_i)
#          = (S_i + T*S_i + S_i+T*S_i) / (S_i / R_insert) = 2 * (1+T) * R_insert
io_rate_per_merge_bps = 2 * (1 + size_ratio_T) * insert_rate_bps

# Total I/O rate is the sum of the initial flush and all subsequent disk merges
total_io_rate_bps = io_rate_mem_flush_bps + (num_disk_to_disk_merges * io_rate_per_merge_bps)

# --- Step 3: Calculate Final Page I/O Rate ---

# Convert the total I/O rate from bytes/sec to pages/sec
total_page_io_rate = total_io_rate_bps / page_size_bytes

print("Calculation of Total I/O Rate in pages/s:")
# The final equation prints all the intermediate numbers used in the calculation
final_equation = f"({io_rate_mem_flush_bps} + {num_disk_to_disk_merges} * 2 * (1 + {int(size_ratio_T)}) * {insert_rate_bps}) / {page_size_bytes}"
print(f"Minimum total page I/O rate = {final_equation}")
print(f"                             = {total_page_io_rate} pages/s")

# <<<876.8>>>