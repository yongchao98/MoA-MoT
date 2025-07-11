# Description: This script calculates the memory consumption per voxel in an optimized smoke simulation.

# --- Initial Full Precision Scheme Analysis ---
# Total memory per voxel in the full precision scheme.
total_bytes_full = 84

# Velocity is twelve 32-bit floating-point numbers. 32 bits = 4 bytes.
num_velocity_vars = 12
bytes_per_float32 = 4
velocity_bytes_full = num_velocity_vars * bytes_per_float32

# Color is three 32-bit floating-point variables.
num_color_vars = 3
color_bytes_full = num_color_vars * bytes_per_float32

# The remaining bytes are for other data (e.g., density, temperature).
other_bytes_full = total_bytes_full - velocity_bytes_full - color_bytes_full

# --- Optimized Scheme Calculation ---
# For optimization while maintaining sufficient precision:
# - Physical quantities like velocity are reduced to 16-bit half-precision floats (2 bytes).
# - Color is reduced to 8-bit unsigned integers per channel (1 byte).

bytes_per_float16 = 2  # 16-bit half-precision float
bytes_per_uint8 = 1    # 8-bit unsigned integer for color channels

# Calculate optimized memory for velocity.
velocity_bytes_optimized = num_velocity_vars * bytes_per_float16

# Calculate optimized memory for color.
color_bytes_optimized = num_color_vars * bytes_per_uint8

# Assume 'other' data fields are also reduced to 16-bit floats.
# First, find how many 'other' variables there were.
num_other_vars = other_bytes_full // bytes_per_float32
# Then, calculate their optimized size.
other_bytes_optimized = num_other_vars * bytes_per_float16

# --- Final Calculation ---
# Sum the components to find the total memory consumption per voxel in the optimized scheme.
total_bytes_optimized = velocity_bytes_optimized + color_bytes_optimized + other_bytes_optimized

# Print the final calculation, showing each optimized component.
# The components are: Optimized Velocity + Optimized Color + Optimized Other Data
print(f"The optimized memory consumption per voxel is the sum of memory for velocity, color, and other properties.")
print(f"Final Calculation: {velocity_bytes_optimized} bytes (Velocity) + {color_bytes_optimized} bytes (Color) + {other_bytes_optimized} bytes (Other) = {total_bytes_optimized} bytes")