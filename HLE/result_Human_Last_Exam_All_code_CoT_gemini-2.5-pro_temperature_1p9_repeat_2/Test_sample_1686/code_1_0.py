# Plan:
# 1. Define the parameters of the full precision scheme.
# 2. Deconstruct the total bytes to find the size of 'other' data.
# 3. Define the optimized data sizes based on standard practices (float32->float16, float32->uint8).
# 4. Calculate the new total memory consumption per voxel.
# 5. Print the final calculation and result.

# --- Step 1 & 2: Full Precision Analysis ---

# A 32-bit float is 4 bytes.
bytes_per_float32 = 4

# Velocity data
num_velocity_vars = 12
full_precision_velocity_bytes = num_velocity_vars * bytes_per_float32

# Color data
num_color_vars = 3
full_precision_color_bytes = num_color_vars * bytes_per_float32

# Total initial data
total_full_precision_bytes = 84

# Calculate the size of the remaining 'other' data
full_precision_other_bytes = total_full_precision_bytes - full_precision_velocity_bytes - full_precision_color_bytes
num_other_vars = full_precision_other_bytes // bytes_per_float32

# --- Step 3: Optimization Scheme ---

# A 16-bit half-precision float is 2 bytes.
bytes_per_float16 = 2
# An 8-bit unsigned integer is 1 byte.
bytes_per_uint8 = 1

# Optimize velocity from float32 to float16
optimized_velocity_bytes = num_velocity_vars * bytes_per_float16

# Optimize color from float32 to uint8 per channel
optimized_color_bytes = num_color_vars * bytes_per_uint8

# Optimize 'other' data from float32 to float16
optimized_other_bytes = num_other_vars * bytes_per_float16


# --- Step 4 & 5: Final Calculation ---

# Sum the components of the optimized scheme
total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

# Print the final equation showing how the total optimized memory is calculated
print("The resulting memory consumption per voxel can be calculated as follows:")
print(f"Optimized Memory = (Velocity) + (Color) + (Other Data)")
print(f"{total_optimized_bytes} bytes = {optimized_velocity_bytes} bytes + {optimized_color_bytes} bytes + {optimized_other_bytes} bytes")
