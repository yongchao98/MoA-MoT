# Initial state analysis
bytes_per_32_bit_float = 4
num_velocity_vars = 12
num_color_vars = 3
total_initial_bytes = 84

# Calculate memory for known initial components
initial_velocity_bytes = num_velocity_vars * bytes_per_32_bit_float
initial_color_bytes = num_color_vars * bytes_per_32_bit_float

# Infer memory and count of other properties
initial_other_bytes = total_initial_bytes - initial_velocity_bytes - initial_color_bytes
num_other_vars = initial_other_bytes // bytes_per_32_bit_float

# Optimization choices (byte sizes)
bytes_per_16_bit_half = 2
bytes_per_8_bit_uint = 1

# Calculate optimized memory consumption
optimized_velocity_bytes = num_velocity_vars * bytes_per_16_bit_half
optimized_color_bytes = num_color_vars * bytes_per_8_bit_uint
optimized_other_bytes = num_other_vars * bytes_per_16_bit_half

# Calculate the final total
total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

# Print the final equation and the result
print("Optimized Memory Calculation:")
print(f"Velocity ({num_velocity_vars} half-floats) + Color ({num_color_vars} bytes) + Other ({num_other_vars} half-floats)")
print(f"Equation: {optimized_velocity_bytes} bytes + {optimized_color_bytes} bytes + {optimized_other_bytes} bytes = {total_optimized_bytes} bytes")
print(f"\nThe resulting memory consumption per voxel would be {total_optimized_bytes} bytes.")

# Final answer in the specified format
print(f"\n<<<{total_optimized_bytes}>>>")