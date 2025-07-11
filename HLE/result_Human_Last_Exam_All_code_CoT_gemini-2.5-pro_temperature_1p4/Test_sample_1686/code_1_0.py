# Number of components for each attribute in a voxel
num_velocity_components = 12
num_color_components = 3
# Derived from the initial state: (84 total - 48 velocity - 12 color) / 4 bytes_per_float = 6
num_other_components = 6

# Optimized byte size for each data type
# Using 16-bit half-precision float for velocity and other fields
bytes_per_velocity_component = 2
bytes_per_other_component = 2
# Using 8-bit unsigned integer for each color channel
bytes_per_color_component = 1

# Calculate the memory consumption for each part of the optimized voxel
velocity_memory = num_velocity_components * bytes_per_velocity_component
color_memory = num_color_components * bytes_per_color_component
other_memory = num_other_components * bytes_per_other_component

# Calculate the total memory consumption
total_memory = velocity_memory + color_memory + other_memory

# Print the final equation and the result
print("Optimized memory consumption calculation:")
print(f"({num_velocity_components} velocity × {bytes_per_velocity_component} bytes) + ({num_color_components} color × {bytes_per_color_component} byte) + ({num_other_components} other × {bytes_per_other_component} bytes) = {total_memory} bytes")