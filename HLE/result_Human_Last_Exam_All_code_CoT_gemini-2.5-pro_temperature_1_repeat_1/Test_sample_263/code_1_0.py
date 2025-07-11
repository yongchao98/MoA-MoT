# The size of the write buffer in kilobytes (KB).
write_buffer_size_kb = 16

# The number of entries that fit in the write buffer.
# The problem statement "the number of entries is 4096" is interpreted
# as the capacity of the write buffer.
num_entries = 4096

# Conversion factor from KB to bytes.
bytes_per_kb = 1024

# Step 1: Calculate the total size of the write buffer in bytes.
write_buffer_size_bytes = write_buffer_size_kb * bytes_per_kb

# Step 2: Calculate the size of a single entry by dividing the total buffer
# size in bytes by the number of entries it can hold.
entry_size_bytes = write_buffer_size_bytes / num_entries

# Final step: Output the equation with the calculated values.
# The final equation shows the division of the buffer size in bytes by the number of entries.
print(f"{write_buffer_size_bytes} / {num_entries} = {int(entry_size_bytes)}")