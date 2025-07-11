import math

# Step 1: Define the given parameters.
# "LSM tree with 6 levels" implies 6 levels on disk (L1 to L6) plus one memory buffer (L0).
num_disk_levels = 6  # L
# Largest level size is 1GB. Using decimal powers for simplicity (10^9).
size_largest_level_bytes = 1 * (10**9)
# Memory buffer is 1KB. Using decimal powers (10^3).
size_buffer_bytes = 1 * (10**3)
# Insert rate in bytes per second.
insert_rate_bytes_per_sec = 16000  # R_insert
# Page size in bytes.
page_size_bytes = 2500  # P

# Step 2: Calculate the size ratio (T).
# The relationship is Size_L = T^L * Size_0.
# T^L = size_largest_level_bytes / size_buffer_bytes
size_ratio_pow_l = size_largest_level_bytes / size_buffer_bytes
# T = (size_ratio_pow_l)^(1/L)
size_ratio = size_ratio_pow_l**(1 / num_disk_levels)
# Round to the nearest integer as T is typically a whole number.
T = int(round(size_ratio))

# Step 3 & 4: Calculate the total I/O rate in bytes per second using the leveled compaction model.
# The total I/O amplification factor for leveled compaction is (2*T*L - 2*T + L).
# This accounts for the initial flush (write) and all subsequent compaction I/O (reads and writes).
amplification_factor = (2 * T * num_disk_levels) - (2 * T) + num_disk_levels
total_io_rate_bytes_per_sec = insert_rate_bytes_per_sec * amplification_factor

# Step 5: Convert the byte rate to page rate.
total_page_io_rate = total_io_rate_bytes_per_sec / page_size_bytes

# Print the final equation and the result.
# The formula for total I/O rate is R_insert * (2*T*L - 2*T + L)
# The final page I/O rate is this value divided by the page size.
print("Problem Parameters:")
print(f"  Number of disk levels (L): {num_disk_levels}")
print(f"  Size Ratio (T): {T}")
print(f"  Insert Rate (R_insert): {insert_rate_bytes_per_sec} bytes/s")
print(f"  Page Size (P): {page_size_bytes} bytes")
print("\nCalculation:")
# The problem asks to output each number in the final equation.
final_equation = f"({insert_rate_bytes_per_sec} * (2 * {T} * {num_disk_levels} - 2 * {T} + {num_disk_levels})) / {page_size_bytes}"
print(f"Minimum total page I/O rate = {final_equation}")
print(f"Result = {total_page_io_rate:.1f} pages/s")
