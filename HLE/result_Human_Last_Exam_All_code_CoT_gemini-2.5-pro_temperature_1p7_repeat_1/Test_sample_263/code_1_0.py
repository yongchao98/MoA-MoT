import math

# Given parameters
levels = 5
size_ratio = 3
total_entries = 4096
buffer_size_kb = 16

# Convert buffer size from KB to bytes
buffer_size_bytes = buffer_size_kb * 1024

# Calculate the multiplier for the geometric series of level sizes
# Multiplier = (T^L - 1) / (T - 1)
series_sum_multiplier = (size_ratio**levels - 1) // (size_ratio - 1)

# The total number of entries in the tree is N_total = N_0 * multiplier,
# where N_0 is the number of entries in the buffer.
# N_0 = S_0_bytes / E, where E is the entry size.
# So, N_total = (S_0_bytes / E) * multiplier
# Solving for E: E = (S_0_bytes * multiplier) / N_total
entry_size = (buffer_size_bytes * series_sum_multiplier) / total_entries

# Output the final equation with all numbers and the result
print("The minimum size of an entry is calculated by the formula:")
print("Entry Size = (Buffer Size in Bytes * Sum of Ratios) / Total Entries")
print("Where Sum of Ratios = (Size Ratio ^ Levels - 1) / (Size Ratio - 1)")
print(f"\n({buffer_size_bytes} * {series_sum_multiplier}) / {total_entries} = {int(entry_size)}")