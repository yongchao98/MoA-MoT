# Initial memory configuration
initial_total_bytes = 84
velocity_float_count = 12
color_float_count = 3
bytes_per_32bit_float = 4

# Calculate the size of each component in the initial scheme
initial_velocity_bytes = velocity_float_count * bytes_per_32bit_float
initial_color_bytes = color_float_count * bytes_per_32bit_float
initial_other_data_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes

# Define the sizes for the optimized data types
bytes_per_16bit_half = 2
bytes_per_8bit_integer = 1

# Calculate the new optimized size for each component
# Velocity and other data are reduced to 16-bit half-precision floats.
optimized_velocity_bytes = velocity_float_count * bytes_per_16bit_half
optimized_other_data_bytes = initial_other_data_bytes / 2

# Color is reduced to 8-bit integers per channel.
optimized_color_bytes = color_float_count * bytes_per_8bit_integer

# Calculate the total memory consumption after optimization
final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_data_bytes

# Print the breakdown of the final calculation
print(f"The resulting memory consumption per voxel would be calculated as follows:")
print(f"Optimized Velocity (16-bit) + Optimized Color (8-bit) + Optimized Other Data (16-bit)")
print(f"{int(optimized_velocity_bytes)} bytes + {int(optimized_color_bytes)} bytes + {int(optimized_other_data_bytes)} bytes = {int(final_total_bytes)} bytes")
