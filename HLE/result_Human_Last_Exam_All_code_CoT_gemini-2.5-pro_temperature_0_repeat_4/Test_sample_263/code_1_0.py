# Define the given parameters for the LSM tree
levels = 5
size_ratio = 3
total_entries = 4096
write_buffer_kb = 16

# Convert the write buffer size from kilobytes to bytes
write_buffer_bytes = write_buffer_kb * 1024

# The formula for the minimum entry size is derived from the total capacity formula of an LSM tree.
# Entry Size = (Buffer_Size_Bytes * (T^L - 1)) / (Total_Entries * (T - 1))
# where L is levels, T is size_ratio.

# Calculate the numerator of the equation
# (Buffer_Size_Bytes * (T^L - 1))
numerator_val = write_buffer_bytes * (size_ratio**levels - 1)

# Calculate the denominator of the equation
# (Total_Entries * (T - 1))
denominator_val = total_entries * (size_ratio - 1)

# Calculate the final entry size
entry_size = numerator_val / denominator_val

# Print the calculation steps with the actual numbers
print("The equation to find the minimum entry size is:")
print(f"Minimum Entry Size = ({write_buffer_bytes} * ({size_ratio}^{levels} - 1)) / ({total_entries} * ({size_ratio} - 1))")
print(f"Minimum Entry Size = ({write_buffer_bytes} * {size_ratio**levels - 1}) / ({total_entries} * {size_ratio - 1})")
print(f"Minimum Entry Size = {int(numerator_val)} / {int(denominator_val)}")
print(f"Minimum Entry Size = {int(entry_size)}")