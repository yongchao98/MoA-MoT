import math

# Step 1: Define the problem parameters
# Number of levels in the LSM tree
num_levels = 6
# Largest level size in bytes (1 GB = 1024^3 bytes)
max_level_size_bytes = 1 * 1024 * 1024 * 1024
# Memory buffer size in bytes (1 KB = 1024 bytes)
mem_buffer_size_bytes = 1 * 1024
# Insert rate in bytes per second
insert_rate_bytes_per_sec = 16000
# Page size in bytes
page_size_bytes = 2500

# Step 2: Calculate the size ratio (T)
# Formula: MaxLevelSize = T^L * MemBufferSize  =>  T = (MaxLevelSize / MemBufferSize)^(1/L)
size_ratio_t = (max_level_size_bytes / mem_buffer_size_bytes)**(1/num_levels)

# Step 3: Calculate the total I/O rate in bytes per second
# I/O amplification factor = L (writes) + (L-1) (reads of data) + (L-1)*T (reads for merging)
# Amplification = L + (L-1)*(1+T)
io_amplification_factor = num_levels + (num_levels - 1) * (1 + size_ratio_t)
total_io_rate_bytes_per_sec = insert_rate_bytes_per_sec * io_amplification_factor

# Step 4: Calculate the minimum total page I/O rate
page_io_rate = total_io_rate_bytes_per_sec / page_size_bytes

# Print the final equation with all the numbers and the result
print(f"The minimum total page I/O rate is calculated as follows:\n")
print(f"Page I/O Rate = (Insert Rate * I/O Amplification Factor) / Page Size")
print(f"I/O Amplification Factor = L + (L - 1) * (1 + T)")
print(f"where L = {num_levels}, T = {size_ratio_t:.4f}")
print("\nFinal equation with values plugged in:")
print(f"Page I/O Rate = ({insert_rate_bytes_per_sec} * ({num_levels} + ({num_levels - 1}) * (1 + {size_ratio_t:.4f}))) / {page_size_bytes}")

print(f"\nResult: {page_io_rate:.3f} pages/s")
