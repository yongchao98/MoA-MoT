import math

# Step 1: Define the given parameters
num_levels_total = 6
mem_buffer_size_kb = 1
largest_level_size_gb = 1
insert_rate_bps = 16000  # bytes per second
page_size_bytes = 2500

# Step 2: Determine LSM tree structure and calculate the size ratio (T)
# We interpret "6 levels" as 1 memory buffer (L0) and 5 disk levels (L1 to L5).
num_disk_levels = num_levels_total - 1

# Convert sizes to bytes
mem_buffer_size_bytes = mem_buffer_size_kb * 1024
largest_level_size_bytes = largest_level_size_gb * (1024**3)

# The size of the largest level L_k is T^k times the size of the initial level L0.
# largest_level_size = T^(num_disk_levels) * mem_buffer_size
# So, T = (largest_level_size / mem_buffer_size)^(1 / num_disk_levels)
size_ratio_base = largest_level_size_bytes / mem_buffer_size_bytes
T = round(math.pow(size_ratio_base, 1 / num_disk_levels))

print("### Step 1: Calculate Size Ratio (T) ###")
print(f"Number of total levels: {num_levels_total}")
print(f"Number of disk levels (L_disk): {num_disk_levels}")
print(f"Memory buffer size (S0): {mem_buffer_size_bytes} bytes")
print(f"Largest level size (S{num_disk_levels}): {largest_level_size_bytes} bytes")
print(f"Size Ratio T = (S{num_disk_levels} / S0) ^ (1/L_disk)")
print(f"T = ({largest_level_size_bytes} / {mem_buffer_size_bytes}) ^ (1/{num_disk_levels}) = {size_ratio_base} ^ {1/num_disk_levels} = {T}")
print("-" * 30)

# Step 3: Calculate the total I/O rate in bytes per second
# Total I/O = (Flush I/O) + (Compaction I/O)
# Flush I/O Rate (L0 -> L1) = insert_rate (writes only)
# Compaction I/O Rate for one merge (Li -> Li+1) = insert_rate * (1_read_Li + T_reads_Li+1 + (T+1)_writes_Li+1)
# Compaction I/O = insert_rate * (1 + T + T + 1) = insert_rate * 2 * (T + 1)
# This occurs for (num_disk_levels - 1) merges.
num_compactions = num_disk_levels - 1
flush_io_rate = insert_rate_bps
compaction_io_rate = num_compactions * insert_rate_bps * 2 * (T + 1)
total_io_rate_bps = flush_io_rate + compaction_io_rate

print("### Step 2: Calculate Total I/O Rate in bytes/s ###")
print("Formula: Total I/O = Flush_I/O + Compaction_I/O")
print("Flush_I/O (L0->L1) = insert_rate")
print("Compaction_I/O (all other levels) = (L_disk - 1) * insert_rate * 2 * (T + 1)")
print(f"Total I/O Rate (bytes/s) = {insert_rate_bps} + ({num_disk_levels} - 1) * {insert_rate_bps} * 2 * ({T} + 1)")
calc_part_1 = num_disk_levels - 1
calc_part_2 = T + 1
calc_part_3 = calc_part_1 * insert_rate_bps * 2 * calc_part_2
print(f"Total I/O Rate (bytes/s) = {insert_rate_bps} + {calc_part_1} * {insert_rate_bps} * 2 * {calc_part_2}")
print(f"Total I/O Rate (bytes/s) = {insert_rate_bps} + {calc_part_3}")
print(f"Total I/O Rate (bytes/s) = {total_io_rate_bps}")
print("-" * 30)

# Step 4: Convert to page I/O rate
total_page_io_rate = total_io_rate_bps / page_size_bytes

print("### Step 3: Calculate Total Page I/O Rate ###")
print("Formula: Page I/O Rate = Total I/O Rate (bytes/s) / Page Size (bytes)")
print(f"Total Page I/O Rate = {total_io_rate_bps} / {page_size_bytes}")
print(f"Final Answer = {total_page_io_rate} pages/s")
