import math

# --- Given Parameters ---
num_levels = 5
size_ratio = 3
total_entries = 4096
buffer_size_kb = 16

# --- Calculations ---

# Convert buffer size from KB to bytes
buffer_size_bytes = buffer_size_kb * 1024

# The total number of entries (N_total) in an L-level LSM tree with size ratio T is
# the sum of entries in all levels. If N_0 is the number of entries in the buffer (level 0),
# then the number of entries in level i is N_0 * T^i.
# N_total = N_0 * (T^0 + T^1 + ... + T^(L-1))
# This is a geometric series with sum: N_total = N_0 * (T^L - 1) / (T - 1)

# Calculate the term (T^L - 1) / (T - 1)
geometric_series_sum = (size_ratio**num_levels - 1) / (size_ratio - 1)

# We can express the number of entries in the buffer (N_0) as:
# N_0 = N_total / geometric_series_sum

# The size of the buffer is also given by:
# Buffer_Size_Bytes = N_0 * Entry_Size_Bytes

# By substituting the expression for N_0, we can solve for the entry size:
# Entry_Size_Bytes = Buffer_Size_Bytes / N_0
# Entry_Size_Bytes = Buffer_Size_Bytes / (N_total / geometric_series_sum)
# Entry_Size_Bytes = (Buffer_Size_Bytes * geometric_series_sum) / N_total

entry_size_bytes = (buffer_size_bytes * geometric_series_sum) / total_entries

# --- Output the Explanation and Result ---

print("Step 1: Determine the total capacity of the LSM tree relative to its write buffer.")
print(f"The ratio of total entries to buffer entries is given by the geometric series sum (T^L - 1) / (T - 1).")
print(f"Calculation: ({size_ratio}^{num_levels} - 1) / ({size_ratio} - 1) = ({int(size_ratio**num_levels)} - 1) / {int(size_ratio - 1)} = {int(geometric_series_sum)}")
print("-" * 20)

print("Step 2: Relate the total number of entries to the number of entries in the write buffer (N_0).")
print(f"Total Entries = N_0 * {int(geometric_series_sum)}")
print(f"{total_entries} = N_0 * {int(geometric_series_sum)}")
print("-" * 20)

print("Step 3: Relate the write buffer size to the size of a single entry.")
print(f"The write buffer size is {buffer_size_kb} KB, which is {buffer_size_kb} * 1024 = {buffer_size_bytes} bytes.")
print(f"Buffer Size = N_0 * Entry Size")
print(f"{buffer_size_bytes} = N_0 * Entry Size")
print("-" * 20)

print("Step 4: Solve for the entry size by combining the equations.")
print("Entry Size = Buffer Size / N_0")
print(f"Substituting for N_0: Entry Size = Buffer Size / (Total Entries / {int(geometric_series_sum)})")
print("This simplifies to the final equation:")
print(f"Entry Size = (Buffer Size * {int(geometric_series_sum)}) / Total Entries")
print(f"Entry Size = ({buffer_size_bytes} * {int(geometric_series_sum)}) / {total_entries}")

print("\nFinal Result:")
print(f"The minimum size of an entry is {int(entry_size_bytes)} bytes.")
