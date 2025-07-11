# Step 1: Define initial data structure based on the problem description.
total_bytes_initial = 84
bytes_per_float32 = 4  # 32-bit float

# Initial components explicitly mentioned
num_velocity_vars = 12
num_color_vars = 3

# Calculate initial size of known components
velocity_bytes_initial = num_velocity_vars * bytes_per_float32
color_bytes_initial = num_color_vars * bytes_per_float32

# Step 2: Determine the size of the other unspecified data.
# The problem states the total is 84 bytes.
other_data_bytes_initial = total_bytes_initial - velocity_bytes_initial - color_bytes_initial
# Assuming this other data is also stored as 32-bit floats
num_other_vars = int(other_data_bytes_initial / bytes_per_float32)

# Step 3: Define optimized data sizes based on a standard optimization scheme.
bytes_per_float16 = 2  # 16-bit half-precision float for simulation data
bytes_per_uint8 = 1    # 8-bit integer for color data

# Step 4: Calculate the new memory consumption for each component.
# Velocity is optimized to half-precision floats.
velocity_bytes_optimized = num_velocity_vars * bytes_per_float16
# Color is optimized to 8-bit integers.
color_bytes_optimized = num_color_vars * bytes_per_uint8
# Other simulation data is also optimized to half-precision floats.
other_data_bytes_optimized = num_other_vars * bytes_per_float16

# Step 5: Calculate the final total size and print the result.
total_bytes_optimized = velocity_bytes_optimized + color_bytes_optimized + other_data_bytes_optimized

print(f"A standard optimization reduces the voxel size by converting simulation data to half-precision (16-bit) floats and color data to 8-bit integers.")
print(f"Optimized Voxel Size Calculation:")
print(f"Velocity   (12 vars * 2 bytes) = {velocity_bytes_optimized} bytes")
print(f"Color      (3 vars * 1 byte)  = {color_bytes_optimized} bytes")
print(f"Other Data (6 vars * 2 bytes)  = {other_data_bytes_optimized} bytes")
print("-" * 35)
# The final required output format, showing each number in the equation
print(f"Final Equation: {velocity_bytes_optimized} + {color_bytes_optimized} + {other_data_bytes_optimized} = {total_bytes_optimized}")
print(f"The resulting memory consumption per voxel is {total_bytes_optimized} bytes.")
