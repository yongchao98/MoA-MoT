# Initial memory layout analysis
initial_total_bytes = 84
bytes_per_float32 = 4

# Velocity component
velocity_float_count = 12
initial_velocity_bytes = velocity_float_count * bytes_per_float32

# Color component
color_float_count = 3
initial_color_bytes = color_float_count * bytes_per_float32

# Unspecified 'other' data (e.g., density, temperature)
initial_other_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
# We can infer the number of float variables used for this 'other' data
other_float_count = initial_other_bytes // bytes_per_float32

print(f"Initial Analysis:")
print(f" - Velocity: {initial_velocity_bytes} bytes ({velocity_float_count} x 32-bit floats)")
print(f" - Color: {initial_color_bytes} bytes ({color_float_count} x 32-bit floats)")
print(f" - Other Data: {initial_other_bytes} bytes ({other_float_count} x 32-bit floats)")
print("-" * 30)

# Optimization assumptions
# - Velocity and Other Data are converted to 16-bit half-precision floats (2 bytes).
# - Color is converted to 8-bit unsigned integers (1 byte per channel).
bytes_per_float16 = 2
bytes_per_uint8 = 1

# Calculate new memory size for each component
optimized_velocity_bytes = velocity_float_count * bytes_per_float16
optimized_color_bytes = color_float_count * bytes_per_uint8
optimized_other_bytes = other_float_count * bytes_per_float16

# Calculate the final total optimized memory size
final_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

print("Optimization Plan:")
print(f" - Velocity will be stored as {velocity_float_count} x 16-bit floats.")
print(f" - Color will be stored as {color_float_count} x 8-bit integers.")
print(f" - Other data will be stored as {other_float_count} x 16-bit floats.")
print("-" * 30)

print("Final Calculation:")
# The request requires printing each number in the final equation.
print(f"Optimized Voxel Size = {optimized_velocity_bytes} (velocity) + {optimized_color_bytes} (color) + {optimized_other_bytes} (other) = {final_optimized_bytes} bytes")
