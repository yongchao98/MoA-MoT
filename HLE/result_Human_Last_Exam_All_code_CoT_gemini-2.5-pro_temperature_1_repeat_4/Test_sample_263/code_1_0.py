import math

# Given parameters
levels = 5
size_ratio = 3
total_entries = 4096
write_buffer_size_kb = 16

# Convert write buffer size from KB to bytes
write_buffer_size_bytes = write_buffer_size_kb * 1024

# The ratio of the total number of entries in the tree to the number of entries
# in the write buffer (level 0) is the sum of a geometric series:
# Ratio = 1 + T + T^2 + ... + T^(L-1) = (T^L - 1) / (T - 1)
total_to_buffer_ratio = (size_ratio**levels - 1) / (size_ratio - 1)

# The size of an entry can be calculated using the following logic:
# 1) total_entries = num_buffer_entries * total_to_buffer_ratio
# 2) write_buffer_size_bytes = num_buffer_entries * entry_size_bytes
# From these, we derive:
# entry_size_bytes = (write_buffer_size_bytes * total_to_buffer_ratio) / total_entries
entry_size_bytes = (write_buffer_size_bytes * total_to_buffer_ratio) / total_entries

print("The formula to calculate the minimum entry size is:")
print("(Write Buffer Size in Bytes * Total-to-Buffer Entry Ratio) / Total Number of Entries")
print("\nPlugging in the numbers:")
print(f"({write_buffer_size_bytes} * {int(total_to_buffer_ratio)}) / {total_entries} = {int(entry_size_bytes)}")
print(f"\nThe minimum size of an entry is {int(entry_size_bytes)} bytes.")