# Plan:
# 1. Define the number of components for velocity, color, and other data based on the problem description.
# 2. Define the byte sizes for the optimized data types:
#    - Half-precision float (16-bit) for velocity and other data = 2 bytes.
#    - Unsigned integer (8-bit) for each color channel = 1 byte.
# 3. Calculate the new memory size in bytes for each component.
# 4. Sum the component sizes to get the final result.
# 5. Print the final equation and the total.

# Number of components from the initial state analysis
num_velocity_components = 12
num_color_components = 3
num_other_data_components = 6 # Deduced from 84 bytes total - (12+3)*4 bytes

# Byte sizes for optimized data types
# Half-precision float (16-bit) for physical properties
half_float_bytes = 2
# Unsigned char/integer (8-bit) for color channels
uchar_bytes = 1

# Calculate the memory for each optimized component
optimized_velocity_bytes = num_velocity_components * half_float_bytes
optimized_color_bytes = num_color_components * uchar_bytes
optimized_other_data_bytes = num_other_data_components * half_float_bytes

# Calculate the total optimized memory consumption
total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_data_bytes

# Print the breakdown of the final calculation and the result
print("Optimized Voxel Memory Calculation:")
print(f"Velocity (12 half-floats)  : {optimized_velocity_bytes} bytes")
print(f"Color (3 unsigned bytes)   : {optimized_color_bytes} bytes")
print(f"Other Data (6 half-floats) : {optimized_other_data_bytes} bytes")
print("---")
print("Total Memory per Voxel:")
# Final equation as requested
print(f"{optimized_velocity_bytes} + {optimized_color_bytes} + {optimized_other_data_bytes} = {total_optimized_bytes} bytes")