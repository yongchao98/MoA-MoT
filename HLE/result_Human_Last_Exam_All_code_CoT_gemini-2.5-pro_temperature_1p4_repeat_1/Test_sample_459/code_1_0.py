import math

# Problem parameters
num_total_levels = 6
largest_level_size_gb = 1
mem_buffer_size_kb = 1
insert_rate_bps = 16000
page_size_bytes = 2500

# Step 1: Define sizes and number of disk levels.
# We interpret the 6 levels as L0 (memory) and L1-L5 (disk).
num_disk_levels = num_total_levels - 1
mem_buffer_size_bytes = mem_buffer_size_kb * 1024
largest_level_size_bytes = largest_level_size_gb * (1024**3)

# Step 2: Calculate the size ratio (T).
# The formula is: largest_level_size = mem_buffer_size * T^(num_disk_levels)
# We solve for T.
T_ratio = (largest_level_size_bytes / mem_buffer_size_bytes)**(1 / num_disk_levels)
# The result should be an integer for a typical LSM tree configuration.
T_ratio = round(T_ratio)

# Step 3: Calculate the total I/O amplification factor for leveled compaction.
# Write amplification: each byte is written once to each disk level.
write_amplification = num_disk_levels
# Read amplification: for each merge (num_disk_levels - 1 merges), data is read from the
# source level (1x) and the destination level (Tx).
read_amplification = (num_disk_levels - 1) * (1 + T_ratio)
# Total amplification is the sum of I/Os per inserted byte.
total_amplification = write_amplification + read_amplification

# Step 4: Calculate the total I/O rate in bytes per second.
total_io_bytes_per_second = total_amplification * insert_rate_bps

# Step 5: Calculate the final page I/O rate.
minimum_page_io_rate = total_io_bytes_per_second / page_size_bytes

# Step 6: Print the final equation with all the numerical values.
# The expanded formula is: ((num_disk_levels + (num_disk_levels - 1) * (1 + T)) * Insert_Rate) / Page_Size
print(f"(({num_disk_levels} + ({num_disk_levels} - 1) * (1 + {int(T_ratio)})) * {insert_rate_bps}) / {page_size_bytes} = {minimum_page_io_rate}")