# LSM Tree parameters
num_levels = 5
size_ratio = 3
num_entries = 4096
write_buffer_size_kb = 16

# Convert write buffer size from KB to bytes
# 1 KB = 1024 bytes
write_buffer_size_bytes = write_buffer_size_kb * 1024

# The size of the write buffer is the product of the number of entries it can hold
# and the size of each entry.
# We can find the size of a single entry by dividing the buffer size by the number of entries.
entry_size_bytes = write_buffer_size_bytes / num_entries

print(f"Given an LSM Tree with:")
print(f"- Number of Entries in Buffer: {num_entries}")
print(f"- Write Buffer Size: {write_buffer_size_kb} KB")
print("\nCalculating the minimum size of an entry...")
print(f"Step 1: Convert buffer size to bytes.")
print(f"   {write_buffer_size_kb} KB * 1024 bytes/KB = {write_buffer_size_bytes} bytes")
print(f"Step 2: Divide buffer size in bytes by the number of entries.")
print(f"   Equation: Entry Size = Total Bytes / Number of Entries")
print(f"   Entry Size = {write_buffer_size_bytes} / {num_entries}")
print(f"\nThe minimum size of an entry is {int(entry_size_bytes)} bytes.")
