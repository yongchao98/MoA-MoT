# 1. Define parameters for the full precision scheme.
float_32_bit_size_bytes = 4
total_initial_bytes = 84

# Number of floating point variables for each property.
velocity_components = 12
color_components = 3

# Calculate initial memory breakdown.
initial_velocity_bytes = velocity_components * float_32_bit_size_bytes
initial_color_bytes = color_components * float_32_bit_size_bytes
initial_other_bytes = total_initial_bytes - initial_velocity_bytes - initial_color_bytes

# Determine the number of other scalar fields.
other_components = initial_other_bytes // float_32_bit_size_bytes

# 2. Define the optimization strategy.
# We will reduce precision from 32-bit floats (4 bytes) to 16-bit half-precision floats (2 bytes).
float_16_bit_size_bytes = 2

# 3. Calculate the new memory consumption for each component.
optimized_velocity_bytes = velocity_components * float_16_bit_size_bytes
optimized_color_bytes = color_components * float_16_bit_size_bytes
optimized_other_bytes = other_components * float_16_bit_size_bytes

# Calculate the total optimized memory.
total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

# Print the final equation with a clear breakdown.
print("Optimized memory calculation per voxel:")
print(f"{optimized_velocity_bytes} bytes (for velocity) + {optimized_color_bytes} bytes (for color) + {optimized_other_bytes} bytes (for other fields) = {total_optimized_bytes} bytes")
