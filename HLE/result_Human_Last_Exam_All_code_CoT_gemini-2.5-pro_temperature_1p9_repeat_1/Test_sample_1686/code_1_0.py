import math

# --- 1. Initial State Analysis ---
# The problem states a total of 84 bytes per voxel.
initial_total_memory_per_voxel = 84

# Breakdown of the known components. A 32-bit float is 4 bytes.
num_velocity_values = 12
num_color_channels = 3
bytes_per_32bit_float = 4

# Calculate memory for specified components.
initial_velocity_mem = num_velocity_values * bytes_per_32bit_float
initial_color_mem = num_color_channels * bytes_per_32bit_float

# The difference between the total and the specified components implies other data is stored.
initial_other_mem = initial_total_memory_per_voxel - initial_velocity_mem - initial_color_mem

# --- 2. Optimization ---
# Propose new, more efficient data types.
# For velocity and other numerical data, use 16-bit half-precision floats (2 bytes).
bytes_per_16bit_float = 2
# For color, use 8-bit unsigned integers (1 byte) per channel.
bytes_per_8bit_uint = 1

# Calculate the new memory footprint for each component.
optimized_velocity_mem = num_velocity_values * bytes_per_16bit_float
optimized_color_mem = num_color_channels * bytes_per_8bit_uint

# Assume the "other" data can also be quantized from 32-bit to 16-bit, halving its memory.
optimized_other_mem = initial_other_mem / 2

# --- 3. Final Calculation ---
# Sum the optimized components to get the final result.
final_total_memory = optimized_velocity_mem + optimized_color_mem + optimized_other_mem

print("An optimized smoke simulation quantization scheme can be achieved as follows:")
print("\n1. Analyze Initial Memory Allocation (per voxel):")
print(f"   - Velocity: {initial_velocity_mem} bytes ({num_velocity_values} values * {bytes_per_32bit_float} bytes)")
print(f"   - Color: {initial_color_mem} bytes ({num_color_channels} channels * {bytes_per_32bit_float} bytes)")
print(f"   - Other Data (Implicit): {int(initial_other_mem)} bytes (derived from total of {initial_total_memory_per_voxel})")

print("\n2. Apply Optimizations:")
print(f"   - Velocity: Convert from 32-bit to 16-bit floats ({bytes_per_32bit_float} -> {bytes_per_16bit_float} bytes/value).")
print(f"   - Color: Convert from 32-bit floats to 8-bit integers ({bytes_per_32bit_float} -> {bytes_per_8bit_uint} byte/channel).")
print(f"   - Other Data: Convert from 32-bit to 16-bit floats (halving memory).")

print("\n3. Calculate Final Memory Consumption (per voxel):")
print("   Memory = (Optimized Velocity) + (Optimized Color) + (Optimized Other Data)")
# Final equation showing each number used in the calculation
print(f"   Final Memory = ({optimized_velocity_mem}) + ({optimized_color_mem}) + ({int(optimized_other_mem)}) = {int(final_total_memory)} bytes")
