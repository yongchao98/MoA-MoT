# Problem parameters
num_levels = 5
size_ratio = 3
total_entries = 4096
buffer_size_kb = 16

# Convert the write buffer size from kilobytes to bytes
buffer_size_bytes = buffer_size_kb * 1024

# The total size of the LSM tree is the sum of a geometric series.
# Formula: S_total = S_0 * (T^L - 1) / (T - 1)
# where S_0 is the size of Level 0 (the buffer)
# T is the size ratio
# L is the number of levels
geometric_sum_factor = (size_ratio**num_levels - 1) / (size_ratio - 1)
total_size_bytes = buffer_size_bytes * geometric_sum_factor

# The minimum entry size is the total size of the tree divided by the total number of entries.
entry_size_bytes = total_size_bytes / total_entries

print(f"To find the minimum entry size, we first calculate the total size of the LSM tree.")
print(f"The total size is the sum of the sizes of all {num_levels} levels.")
print(f"Total Tree Size = (Buffer Size) * (Ratio^Levels - 1) / (Ratio - 1)")
print(f"Total Tree Size = {int(buffer_size_bytes)} bytes * ({size_ratio}^{num_levels} - 1) / ({size_ratio} - 1) = {int(total_size_bytes)} bytes")
print("")
print(f"Next, we divide the total tree size by the total number of entries to find the size per entry.")
print(f"Entry Size = Total Tree Size / Total Number of Entries")
# Final equation with numbers
print(f"Entry Size = {int(total_size_bytes)} / {total_entries} = {int(entry_size_bytes)} bytes")
