import math

# --- Given Parameters ---
levels = 5
size_ratio = 3
total_entries = 4096
buffer_size_kb = 16

# --- Calculations ---

# 1. Convert the write buffer size from kilobytes to bytes
buffer_size_bytes = buffer_size_kb * 1024

# 2. Calculate the total size multiplier
# The multiplier is the sum of the geometric series: 1 + T + T^2 + ... + T^(L-1)
# The formula for this sum is (T^L - 1) / (T - 1)
size_multiplier = (size_ratio**levels - 1) / (size_ratio - 1)

# 3. Calculate the minimum size of an entry
# The equation is: Entry Size = (Buffer Size * Multiplier) / Total Entries
entry_size = (buffer_size_bytes * size_multiplier) / total_entries

# --- Output the result ---
print(f"To find the minimum size of an entry, we set up the following equation:")
print("Entry Size = (Write Buffer Size in Bytes * Size Multiplier) / Total Entries")
print("\nFirst, we calculate the values for the equation:")
print(f"Write Buffer Size = {buffer_size_kb}KB = {buffer_size_bytes} bytes")
print(f"Size Multiplier = ({size_ratio}^{levels} - 1) / ({size_ratio} - 1) = {int(size_multiplier)}")
print(f"Total Entries = {total_entries}")

print("\nPlugging these numbers into the final equation, we get:")
# Using math.trunc to ensure the output is an integer for display
print(f"{math.trunc(entry_size)} = ({buffer_size_bytes} * {math.trunc(size_multiplier)}) / {total_entries}")
print(f"\nThe minimum size of an entry is {math.trunc(entry_size)} bytes.")
