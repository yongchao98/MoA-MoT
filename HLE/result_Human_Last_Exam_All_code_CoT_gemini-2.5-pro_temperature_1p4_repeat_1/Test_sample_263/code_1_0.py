import math

# --- Given Parameters ---
num_levels = 5
size_ratio = 3
total_entries = 4096
buffer_size_kb = 16

# --- Calculation ---

# 1. Convert write buffer size from KB to Bytes.
buffer_size_bytes = buffer_size_kb * 1024

# 2. Calculate the total capacity multiplier for the tree.
# This is derived from the geometric series sum: (T^L - 1) / (T - 1).
capacity_multiplier_numerator = (size_ratio ** num_levels) - 1
capacity_multiplier_denominator = size_ratio - 1
capacity_multiplier = capacity_multiplier_numerator // capacity_multiplier_denominator

# 3. Calculate the total capacity of the tree in bytes.
total_capacity_bytes = buffer_size_bytes * capacity_multiplier

# 4. Calculate the minimum size of a single entry.
min_entry_size = total_capacity_bytes // total_entries

# --- Output ---
# The final output prints the full equation with all the initial values,
# showing how the final result is derived.
print(f"The minimum size of an entry is calculated by the formula:")
print(f"Entry Size = (Buffer Size in Bytes * Capacity Multiplier) / Total Entries")
print(f"\nPlugging in the numbers:")
print(f"(({buffer_size_kb} * 1024) * (({size_ratio} ** {num_levels} - 1) / ({size_ratio} - 1))) / {total_entries} = {min_entry_size}")