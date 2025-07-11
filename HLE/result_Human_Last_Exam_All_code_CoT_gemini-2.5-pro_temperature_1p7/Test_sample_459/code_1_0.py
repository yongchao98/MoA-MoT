import math

# Step 1: Define the given parameters
num_levels = 6
largest_level_size_gb = 1
mem_buffer_size_kb = 1
insert_rate_bytes_per_sec = 16000
page_size_bytes = 2500

# Convert sizes to bytes for consistent units
largest_level_size_bytes = largest_level_size_gb * (1024**3)
mem_buffer_size_bytes = mem_buffer_size_kb * 1024

# Step 2: Calculate the Size Ratio (T)
# The model assumes Size(L_i) = T * Size(L_{i-1}), so LargestLevelSize = T^NumLevels * BufferSize
size_ratio_base = largest_level_size_bytes / mem_buffer_size_bytes
T = size_ratio_base ** (1 / num_levels)

# Step 3: Calculate the total I/O rate in bytes per second
# In a leveled LSM tree, total I/O for each inserted byte is:
# 1 write (memory to L1)
# + (L-1) merges, where each merge costs (1 read_small + T reads_large + 1 write_large)
# Total I/O Multiplier = 1 + (num_levels - 1) * (1 + T + 1)
# Total I/O Multiplier = 1 + (num_levels - 1) * (T + 2)
io_multiplier = 1 + (num_levels - 1) * (T + 2)
total_io_bytes_per_sec = insert_rate_bytes_per_sec * io_multiplier

# Step 4: Calculate the minimum total page I/O rate
total_page_io_rate = total_io_bytes_per_sec / page_size_bytes

# Print the explanation and the final equation
print("--- LSM Tree I/O Calculation ---")
print(f"Given Parameters:")
print(f"  - Number of Disk Levels (L): {num_levels}")
print(f"  - Largest Level Size: {largest_level_size_gb} GB")
print(f"  - Memory Buffer Size: {mem_buffer_size_kb} KB")
print(f"  - Insert Rate: {insert_rate_bytes_per_sec} bytes/s")
print(f"  - Page Size: {page_size_bytes} bytes")
print("-" * 34)

print("Step-by-Step Calculation:")
print(f"1. Size Ratio (T) = ({largest_level_size_bytes:.0f} / {mem_buffer_size_bytes:.0f})^(1/{num_levels}) = {T:.4f}")
print(f"2. Total I/O (bytes/s) = {insert_rate_bytes_per_sec} * (1 + ({num_levels}-1) * ({T:.4f} + 2)) = {total_io_bytes_per_sec:.2f} bytes/s")
print(f"3. Page I/O Rate = {total_io_bytes_per_sec:.2f} / {page_size_bytes} = {total_page_io_rate:.2f} pages/s")
print("-" * 34)

print("Final Equation with calculated values:")
# The equation format is: (Insert Rate * (1 + (L - 1) * (T + 2))) / Page Size
print(f"({insert_rate_bytes_per_sec:.0f} * (1 + ({num_levels} - 1) * ({T:.4f} + 2))) / {page_size_bytes:.0f} = {total_page_io_rate:.2f}")

<<<392.94>>>