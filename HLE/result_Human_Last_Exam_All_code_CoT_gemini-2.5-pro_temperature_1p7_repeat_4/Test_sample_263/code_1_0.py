# Problem parameters
num_levels = 5
size_ratio = 3
total_entries = 4096
write_buffer_size_kb = 16

# 1. Calculate the write buffer size in bytes
bytes_per_kb = 1024
buffer_size_bytes = write_buffer_size_kb * bytes_per_kb

# 2. Calculate the sum factor from the geometric series formula for LSM tree capacity
# Total Entries = N0 * (T^L - 1) / (T - 1)
# Let's call the part after N0 the 'sum_factor'
sum_factor = (size_ratio**num_levels - 1) / (size_ratio - 1)

# 3. Calculate the number of entries in the write buffer (N0)
# N0 = Total Entries / sum_factor
num_entries_in_buffer = total_entries / sum_factor

# 4. Calculate the minimum size of one entry in bytes
# Entry Size = Buffer Size (bytes) / Number of Entries in Buffer
min_entry_size = buffer_size_bytes / num_entries_in_buffer

# Print the final equation with all the numbers
print("The final equation to find the minimum entry size is:")
print("Minimum Entry Size = Write Buffer Size (bytes) / Number of Entries in Buffer")
print("\nSubstituting the values:")
# We show the full calculation in one line to represent the final equation with its components
# Using int() for sum_factor because it's a whole number, to make the output cleaner.
print(f"{buffer_size_bytes} / ({total_entries} / {int(sum_factor)}) = {min_entry_size}")

print(f"\nThe minimum size of an entry is {int(min_entry_size)} bytes.")
