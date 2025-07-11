# Plan:
# 1. Deconstruct the initial memory layout to understand all components.
# 2. Propose a standard optimization by reducing data type precision.
# 3. Calculate the new memory size for each component and sum them up.

# --- Initial State Analysis ---
total_mem_full_precision = 84  # bytes
num_velocity_floats = 12
num_color_floats = 3
bytes_per_32_bit_float = 4

# Calculate memory usage of known components in the initial scheme
mem_velocity_full = num_velocity_floats * bytes_per_32_bit_float
mem_color_full = num_color_floats * bytes_per_32_bit_float

# Calculate the size of the remaining unspecified data
mem_other_data_full = total_mem_full_precision - mem_velocity_full - mem_color_full

# --- Optimization ---
# Propose new data types that save memory while maintaining sufficient precision.
# For velocity and other data: Use 16-bit half-precision floats (2 bytes).
# For color: Use 8-bit unsigned integers (1 byte) for each of the 3 channels.
bytes_per_16_bit_float = 2
bytes_per_8_bit_uint = 1

# Calculate memory usage of each component after optimization
mem_velocity_optimized = num_velocity_floats * bytes_per_16_bit_float
# The 3 color channels will now each use 1 byte
mem_color_optimized = 3 * bytes_per_8_bit_uint
# Assume the "other data" is also float-based and can be halved in size
mem_other_data_optimized = mem_other_data_full / 2

# Calculate the total memory consumption after optimization
total_mem_optimized = mem_velocity_optimized + mem_color_optimized + mem_other_data_optimized

# --- Output the final result ---
print("Optimized Memory Calculation per Voxel")
print("--------------------------------------")
print(f"Optimized Velocity Memory (12 * 16-bit floats): {mem_velocity_optimized} bytes")
print(f"Optimized Color Memory (3 * 8-bit integers): {mem_color_optimized} bytes")
print(f"Optimized Other Data Memory (halved from {mem_other_data_full} bytes): {int(mem_other_data_optimized)} bytes")
print("--------------------------------------")
# Final equation as requested, showing each number
print(f"Total Optimized Memory = {mem_velocity_optimized} bytes (velocity) + {mem_color_optimized} bytes (color) + {int(mem_other_data_optimized)} bytes (other) = {int(total_mem_optimized)} bytes")