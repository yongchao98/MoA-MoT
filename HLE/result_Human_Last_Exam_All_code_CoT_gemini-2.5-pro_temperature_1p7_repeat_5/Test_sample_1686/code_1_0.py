# Initial memory configuration analysis
total_initial_bytes = 84
bytes_per_32_bit_float = 32 // 8

# Velocity component
num_velocity_floats = 12
velocity_initial_bytes = num_velocity_floats * bytes_per_32_bit_float

# Color component
num_color_channels = 3
color_initial_bytes = num_color_channels * bytes_per_32_bit_float

# Unspecified data component
unspecified_initial_bytes = total_initial_bytes - velocity_initial_bytes - color_initial_bytes
num_unspecified_floats = unspecified_initial_bytes // bytes_per_32_bit_float

# Optimized memory configuration based on standard quantization schemes
bytes_per_16_bit_float = 16 // 8
bytes_per_8_bit_integer = 8 // 8

# Calculate memory for optimized components
# Velocity and other physical data are reduced to half-precision (16-bit) floats
velocity_optimized_bytes = num_velocity_floats * bytes_per_16_bit_float
unspecified_optimized_bytes = num_unspecified_floats * bytes_per_16_bit_float

# Color data is quantized to 8-bit integers per channel
color_optimized_bytes = num_color_channels * bytes_per_8_bit_integer

# Calculate total optimized memory per voxel
total_optimized_bytes = velocity_optimized_bytes + color_optimized_bytes + unspecified_optimized_bytes

# Print the final equation showing the size of each component
print("Optimized memory per voxel is the sum of memory for each optimized component:")
print(f"Velocity (bytes) + Color (bytes) + Unspecified (bytes) = Total (bytes)")
print(f"{int(velocity_optimized_bytes)} + {int(color_optimized_bytes)} + {int(unspecified_optimized_bytes)} = {int(total_optimized_bytes)}")