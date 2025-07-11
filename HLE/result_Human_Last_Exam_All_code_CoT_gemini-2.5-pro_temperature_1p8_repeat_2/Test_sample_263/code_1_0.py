# Define the given parameters
L = 5  # number of levels
T = 3  # size ratio
N = 4096  # total number of entries
write_buffer_size_kb = 16 # write buffer size in KB

# Convert write buffer size to bytes
write_buffer_size_bytes = write_buffer_size_kb * 1024

# The formula for the minimum size of an entry is derived from the structure of the LSM tree.
# Entry_Size = (Write_Buffer_Size_in_bytes * (T^L - 1)) / (N_total_entries * (T - 1))

# Calculate the intermediate components of the formula
t_pow_l_minus_1 = T**L - 1
t_minus_1 = T - 1

# Calculate the numerator and denominator of the final equation
numerator = write_buffer_size_bytes * t_pow_l_minus_1
denominator = N * t_minus_1

# Calculate the final result
entry_size = numerator / denominator

# Output the step-by-step calculation with the numbers plugged in
print("To find the minimum entry size, we use the following formula:")
print("Entry_Size = (Write_Buffer_Size_Bytes * (T^L - 1)) / (N_total * (T - 1))")
print("\nPlugging in the given values:")
print(f"Entry_Size = ({write_buffer_size_bytes} * ({T}^{L} - 1)) / ({N} * ({T} - 1))")
print(f"Entry_Size = ({write_buffer_size_bytes} * {t_pow_l_minus_1}) / ({N} * {t_minus_1})")
print(f"Entry_Size = {int(numerator)} / {int(denominator)}")
print(f"\nThe minimum size of an entry is {int(entry_size)} bytes.")

print(f"\n<<<{int(entry_size)}>>>")