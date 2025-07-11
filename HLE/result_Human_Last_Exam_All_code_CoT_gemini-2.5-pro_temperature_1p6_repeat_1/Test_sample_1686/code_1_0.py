# Step 1: Define parameters for the initial full precision scheme.
total_bytes_initial = 84
bytes_per_float32 = 4  # 32-bit float

# Calculate initial memory for specified components.
num_velocity_values = 12
velocity_bytes_initial = num_velocity_values * bytes_per_float32

num_color_values = 3
color_bytes_initial = num_color_values * bytes_per_float32

# Calculate the size of the unspecified "other data".
other_bytes_initial = total_bytes_initial - velocity_bytes_initial - color_bytes_initial

# Step 2: Define parameters for the proposed optimized scheme.
bytes_per_float16 = 2  # 16-bit half-precision float
bytes_per_uint8 = 1    # 8-bit unsigned integer for color channels

# Step 3: Calculate the memory for each component in the optimized scheme.

# Optimize velocity from 32-bit float to 16-bit float.
velocity_bytes_optimized = num_velocity_values * bytes_per_float16

# Optimize color from 32-bit float to 8-bit integer.
color_bytes_optimized = num_color_values * bytes_per_uint8

# Assume "other data" was also stored as 32-bit floats and optimize to 16-bit floats.
# First, determine how many values were in "other data".
num_other_values = other_bytes_initial // bytes_per_float32
other_bytes_optimized = num_other_values * bytes_per_float16

# Step 4: Calculate the total memory consumption for the optimized scheme.
total_bytes_optimized = velocity_bytes_optimized + color_bytes_optimized + other_bytes_optimized

# Print the breakdown and the final equation.
print("Based on a standard optimization:")
print(f"- Optimized velocity memory: {velocity_bytes_optimized} bytes")
print(f"- Optimized color memory: {color_bytes_optimized} bytes")
print(f"- Optimized other data memory: {other_bytes_optimized} bytes")
print("\nThe final equation for the total memory per voxel is:")
print(f"{velocity_bytes_optimized} + {color_bytes_optimized} + {other_bytes_optimized} = {total_bytes_optimized}")
print(f"\nThe resulting memory consumption per voxel is {total_bytes_optimized} bytes.")
