# Define the given parameters
num_levels = 5      # Number of levels in the LSM tree (contextual info)
size_ratio = 3      # Size ratio between levels (contextual info)
num_entries = 4096  # Number of entries in the write buffer
buffer_size_kb = 16 # Size of the write buffer in KB

# Define the conversion factor from KB to bytes
bytes_per_kb = 1024

# Step 1: Calculate the total size of the write buffer in bytes
buffer_size_bytes = buffer_size_kb * bytes_per_kb

# Step 2: Calculate the size of a single entry in bytes
# This is derived from the formula: Buffer Size = Number of Entries * Entry Size
entry_size_bytes = buffer_size_bytes / num_entries

# Print the explanation and the result
print("To find the minimum size of an entry, we divide the total size of the write buffer by the number of entries it holds.")
print(f"Write Buffer Size: {buffer_size_kb} KB")
print(f"Number of Entries: {num_entries}")
print("")
print(f"First, we convert the buffer size to bytes:")
print(f"{buffer_size_kb} KB * {bytes_per_kb} bytes/KB = {buffer_size_bytes} bytes")
print("")
print("Then, we divide the buffer size in bytes by the number of entries:")
print(f"Final Equation: ({buffer_size_kb} * {bytes_per_kb}) / {num_entries} = {int(entry_size_bytes)}")
print("")
print(f"The minimum size of an entry is {int(entry_size_bytes)} bytes.")
