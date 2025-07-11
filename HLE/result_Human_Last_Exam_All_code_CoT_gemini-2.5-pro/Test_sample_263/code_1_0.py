# LSM Tree Parameters
# Total levels in the tree. This information is not needed for the final calculation.
levels = 5
# Size ratio between levels. This information is not needed for the final calculation.
size_ratio = 3
# Number of entries in the write buffer (Level 0).
num_entries_buffer = 4096
# Size of the write buffer in kilobytes.
buffer_size_kb = 16

# 1. Convert the write buffer size from kilobytes to bytes.
# 1 KB = 1024 bytes.
bytes_per_kb = 1024
buffer_size_bytes = buffer_size_kb * bytes_per_kb

# 2. Calculate the minimum size of a single entry.
# This is the total size of the buffer in bytes divided by the number of entries it holds.
entry_size_bytes = buffer_size_bytes / num_entries_buffer

# 3. Print the final equation showing how the result was calculated.
# The result should be an integer, so we cast it to int for clean printing.
print("The equation to find the minimum entry size is:")
print("Entry Size (bytes) = Buffer Size (bytes) / Number of Entries in Buffer")
print(f"{int(entry_size_bytes)} = {buffer_size_bytes} / {num_entries_buffer}")

# The final answer is the calculated entry size.
# print(f"The minimum size of an entry is {int(entry_size_bytes)} bytes.")