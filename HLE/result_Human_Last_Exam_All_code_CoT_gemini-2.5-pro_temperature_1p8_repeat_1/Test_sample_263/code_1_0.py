# Parameters from the problem statement
num_entries = 4096
buffer_size_kb = 16

# Conversion factor for KB to Bytes
bytes_per_kb = 1024

# Step 1: Calculate the total size of the write buffer in bytes
buffer_size_bytes = buffer_size_kb * bytes_per_kb

# Step 2: Calculate the size of a single entry by dividing the buffer size by the number of entries
# This assumes the "4096 entries" refers to the capacity of the write buffer.
entry_size_bytes = buffer_size_bytes / num_entries

# Step 3: Print the full equation and the final result as requested.
# The final result is cast to an integer, as the calculation results in a whole number.
print(f"To find the size of an entry, we divide the buffer size in bytes by the number of entries.")
print(f"The equation is: ({buffer_size_kb} * {bytes_per_kb}) / {num_entries}")
print(f"Result: {buffer_size_bytes} / {num_entries} = {int(entry_size_bytes)}")
print(f"\nThe minimum size of an entry is {int(entry_size_bytes)} bytes.")
