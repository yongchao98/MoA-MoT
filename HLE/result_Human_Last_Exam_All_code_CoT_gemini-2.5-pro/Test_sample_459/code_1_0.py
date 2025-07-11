# Step 1: Define initial parameters and establish the LSM tree structure.
num_levels_total = 6
size_largest_level_gb = 1.0
size_mem_buffer_kb = 1.0
insert_rate_bps = 16000.0
page_size_bytes = 2500.0

# Assuming 1 memory level and 5 disk levels
num_disk_levels = num_levels_total - 1

# Convert all sizes to a consistent unit (bytes)
size_mem_buffer_bytes = size_mem_buffer_kb * 1024
size_largest_level_bytes = size_largest_level_gb * (1024**3)

print("Step 1: LSM Tree Parameters")
print(f"  - Number of disk levels: {num_disk_levels}")
print(f"  - Memory buffer size: {size_mem_buffer_bytes} bytes")
print(f"  - Largest level size: {size_largest_level_bytes} bytes")
print("-" * 20)

# Step 2: Calculate the size ratio (T)
# Formula: Size(L_k) = T^k * Size_mem, where k is the number of disk levels.
# T = (Size_largest / Size_mem)^(1/k)
size_ratio = (size_largest_level_bytes / size_mem_buffer_bytes)**(1 / num_disk_levels)
# Let's round T to the nearest integer as it's very close to 16.0
T = round(size_ratio)

print("Step 2: Calculate Size Ratio (T)")
print(f"  - T = (Largest Level Size / Memory Buffer Size) ^ (1 / Number of Disk Levels)")
print(f"  - T = ({size_largest_level_bytes} / {size_mem_buffer_bytes}) ^ (1 / {num_disk_levels})")
print(f"  - T = {size_ratio:.4f} ≈ {T}")
print("-" * 20)

# Step 3: Calculate the total I/O rate in bytes per second.
# Total I/O = (Initial Flush) + (I/O from compactions between disk levels)
# Initial Flush Rate = insert_rate_bps (writes only)
# Number of compactions between disk levels = num_disk_levels - 1
# I/O rate per compaction = insert_rate_bps * 2 * (T + 1)
# Total I/O Rate = insert_rate_bps * (1 + (num_disk_levels - 1) * 2 * (T + 1))
num_compactions = num_disk_levels - 1
total_io_rate_bps = insert_rate_bps * (1 + num_compactions * 2 * (T + 1))

print("Step 3: Calculate Total I/O Rate in bytes/s")
print(f"  - Total I/O Rate = Insert Rate * (1 + (Num Compactions) * 2 * (T + 1))")
print(f"  - Total I/O Rate = {insert_rate_bps} * (1 + {num_compactions} * 2 * ({T} + 1))")
print(f"  - Total I/O Rate = {insert_rate_bps} * (1 + {num_compactions * 2 * (T + 1)})")
print(f"  - Total I/O Rate = {total_io_rate_bps:.2f} bytes/s")
print("-" * 20)

# Step 4: Calculate the total page I/O rate.
# Page I/O Rate = Total I/O Rate / Page Size
page_io_rate = total_io_rate_bps / page_size_bytes

print("Step 4: Calculate Total Page I/O Rate")
print(f"  - Page I/O Rate = (Total I/O Rate) / (Page Size)")
print(f"  - Page I/O Rate = {total_io_rate_bps:.2f} / {page_size_bytes}")
print(f"  - Final Answer = {page_io_rate:.2f} pages/s")
print("-" * 20)