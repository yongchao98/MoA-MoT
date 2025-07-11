import math

# Given parameters
levels = 5
size_ratio = 3
num_entries = 4096
buffer_size_kb = 16

# --- Calculation ---

# Step 1: Convert the write buffer size from KB to bytes.
# 1 KB = 1024 bytes.
buffer_size_bytes = buffer_size_kb * 1024

# Step 2: Calculate the total size of the LSM tree using the geometric series sum formula.
# Total Size = S₀ * (Tᴸ - 1) / (T - 1)
# Calculate the intermediate terms for clarity.
term_in_parentheses_numerator = (size_ratio ** levels) - 1
term_in_parentheses_denominator = size_ratio - 1
size_multiplier = term_in_parentheses_numerator / term_in_parentheses_denominator

total_size_bytes = buffer_size_bytes * size_multiplier

# Step 3: Calculate the size of a single entry.
# Entry Size = Total Size / Total Number of Entries
entry_size_bytes = total_size_bytes / num_entries

# --- Output the result and the equation ---

print("To find the minimum size of an entry, we first calculate the total size of the LSM tree and then divide it by the number of entries.")
print("\nHere is the step-by-step calculation:")

print(f"\nFinal Equation:")
print(f"Entry Size = (Write Buffer Size * (Size Ratio ^ Levels - 1) / (Size Ratio - 1)) / Total Entries")
print(f"Entry Size = ({buffer_size_bytes} * ({size_ratio} ^ {levels} - 1) / ({size_ratio} - 1)) / {num_entries}")
print(f"Entry Size = ({buffer_size_bytes} * ({int(math.pow(size_ratio, levels))} - 1) / {term_in_parentheses_denominator}) / {num_entries}")
print(f"Entry Size = ({buffer_size_bytes} * {term_in_parentheses_numerator} / {term_in_parentheses_denominator}) / {num_entries}")
print(f"Entry Size = ({buffer_size_bytes} * {int(size_multiplier)}) / {num_entries}")
print(f"Entry Size = {int(total_size_bytes)} / {num_entries}")
print(f"Entry Size = {int(entry_size_bytes)} bytes")