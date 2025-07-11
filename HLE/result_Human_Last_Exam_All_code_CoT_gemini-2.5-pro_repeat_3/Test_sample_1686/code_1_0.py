import math

# Step 1: Define the initial full precision state.
# A 32-bit float is 4 bytes.
bytes_per_float32 = 4
initial_total_bytes = 84

# Number of variables for known components.
velocity_floats_count = 12
color_floats_count = 3

# Calculate memory for the known components in the full precision scheme.
velocity_bytes_full = velocity_floats_count * bytes_per_float32
color_bytes_full = color_floats_count * bytes_per_float32

# Step 2: Deduce the memory used by other unspecified physical quantities.
other_data_bytes_full = initial_total_bytes - velocity_bytes_full - color_bytes_full
# Determine how many 32-bit float variables this corresponds to.
other_data_floats_count = other_data_bytes_full / bytes_per_float32

# Step 3: Define the optimized precision sizes.
# Physical quantities (velocity, other data) are reduced to 16-bit half-precision floats (2 bytes).
bytes_per_float16 = 2
# Color data is reduced to 8-bit unsigned integers (1 byte per channel).
bytes_per_uint8 = 1

# Step 4: Calculate the new memory size for each component after optimization.
velocity_bytes_optimized = velocity_floats_count * bytes_per_float16
color_bytes_optimized = color_floats_count * bytes_per_uint8
other_data_bytes_optimized = other_data_floats_count * bytes_per_float16

# Calculate the final total memory consumption.
final_total_bytes = velocity_bytes_optimized + color_bytes_optimized + other_data_bytes_optimized

# Print the breakdown and the final equation.
print("Analysis of Memory Consumption per Voxel")
print("----------------------------------------")
print(f"Initial State (Full Precision): {initial_total_bytes} bytes")
print(f" - Velocity: {velocity_bytes_full} bytes ({velocity_floats_count} x 32-bit floats)")
print(f" - Color: {color_bytes_full} bytes ({color_floats_count} x 32-bit floats)")
print(f" - Other Data: {int(other_data_bytes_full)} bytes ({int(other_data_floats_count)} x 32-bit floats)")
print("\nOptimized State (Reduced Precision):")
print(f" - Optimized Velocity: {velocity_bytes_optimized} bytes ({velocity_floats_count} x 16-bit floats)")
print(f" - Optimized Color: {color_bytes_optimized} bytes ({color_floats_count} x 8-bit integers)")
print(f" - Optimized Other Data: {int(other_data_bytes_optimized)} bytes ({int(other_data_floats_count)} x 16-bit floats)")
print("\n----------------------------------------")
print("Final Calculation:")
# Using math.ceil to ensure integer representation in the final printout.
# Although the numbers are clean integers in this case, it's good practice.
final_equation_str = (
    f"{math.ceil(velocity_bytes_optimized)} (velocity) + "
    f"{math.ceil(color_bytes_optimized)} (color) + "
    f"{math.ceil(other_data_bytes_optimized)} (other) = "
    f"{math.ceil(final_total_bytes)} bytes"
)
print(final_equation_str)
print("\nThe resulting memory consumption per voxel is {} bytes.".format(math.ceil(final_total_bytes)))