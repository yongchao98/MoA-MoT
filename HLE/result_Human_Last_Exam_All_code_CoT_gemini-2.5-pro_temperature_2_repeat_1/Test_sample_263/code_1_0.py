import math

# Problem parameters from the user's request
num_levels = 5
size_ratio = 3
total_entries = 4096
write_buffer_kb = 16

# Convert the write buffer size from kilobytes to bytes
write_buffer_bytes = write_buffer_kb * 1024

# --- Calculate the numerator of our derived formula ---
# (Write Buffer Size in Bytes * (Size Ratio ^ Levels - 1))
term_in_numerator_calc = size_ratio**num_levels - 1
final_numerator = write_buffer_bytes * term_in_numerator_calc

# --- Calculate the denominator of our derived formula ---
# (Total Entries * (Size Ratio - 1))
term_in_denominator_calc = size_ratio - 1
final_denominator = total_entries * term_in_denominator_calc

# --- Final Calculation ---
entry_size = final_numerator / final_denominator

# --- Print the explanation and final equation ---
print("The formula for the minimum entry size is derived from the total capacity and buffer size relationships:")
print("Entry Size = (Write Buffer Size in Bytes * (Size Ratio ^ Levels - 1)) / (Total Entries * (Size Ratio - 1))\n")
print("Here is the final equation with all the numbers filled in:\n")

# Use f-string to show the final equation with values plugged in
print(f"Entry Size = ({write_buffer_bytes} * ({size_ratio}^{num_levels} - 1)) / ({total_entries} * ({size_ratio} - 1))")
print(f"Entry Size = ({write_buffer_bytes} * {int(term_in_numerator_calc)}) / ({total_entries} * {term_in_denominator_calc})")
print(f"Entry Size = {int(final_numerator)} / {int(final_denominator)}")

print(f"\nThe minimum size of an entry is {int(entry_size)} bytes.")