import math

# Given parameters for the LSM tree
levels = 5
size_ratio = 3
total_entries = 4096
write_buffer_size_kb = 16

# Convert write buffer size from kilobytes to bytes
write_buffer_size_bytes = write_buffer_size_kb * 1024

# The total number of entries is the sum of entries in all levels.
# This can be expressed as N_total = N_0 * (1 + T + T^2 + ... + T^(L-1)),
# where N_0 is the number of entries in the write buffer (Level 0),
# T is the size ratio, and L is the number of levels.
# The sum is a geometric series: (T^L - 1) / (T - 1).

# Calculate the sum of the geometric series ratio
sum_ratio_numerator = size_ratio**levels - 1
sum_ratio_denominator = size_ratio - 1
sum_ratio = sum_ratio_numerator / sum_ratio_denominator

# From N_total = N_0 * sum_ratio, we can find N_0.
# N_0 = N_total / sum_ratio.

# The size of an entry is the total size of the write buffer divided by the number of entries it holds.
# Entry Size = write_buffer_size_bytes / N_0
# Entry Size = write_buffer_size_bytes / (total_entries / sum_ratio)
# Entry Size = (write_buffer_size_bytes * sum_ratio) / total_entries
entry_size = (write_buffer_size_bytes * sum_ratio) / total_entries

# Print the final equation with all the numbers used in the calculation
print("The equation to find the minimum entry size is:")
print("(Write Buffer Size in Bytes * Sum of Size Ratios) / Total Entries = Entry Size")
print("")
print("Plugging in the numbers:")
# The final equation string showing all the components
final_equation = f"({write_buffer_size_kb} * 1024 * (({size_ratio}^{levels} - 1) / ({size_ratio} - 1))) / {total_entries} = {int(entry_size)}"
print(final_equation)
print("")
# Simplified calculation
print("Simplified equation:")
simplified_equation = f"({write_buffer_size_bytes} * {int(sum_ratio)}) / {total_entries} = {int(entry_size)}"
print(simplified_equation)