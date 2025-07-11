import math

# --- Initial State Analysis ---
# Given total size per voxel
total_initial_bytes = 84

# Size of a single 32-bit float in bytes
float32_size_bytes = 32 / 8

# 1. Calculate initial memory for Velocity
velocity_var_count = 12
velocity_initial_bytes = velocity_var_count * float32_size_bytes

# 2. Calculate initial memory for Color
color_var_count = 3
color_initial_bytes = color_var_count * float32_size_bytes

# 3. Calculate memory for other unspecified data
other_initial_bytes = total_initial_bytes - (velocity_initial_bytes + color_initial_bytes)

# --- Optimization Strategy ---

# 1. Optimize Velocity: Reduce from 32-bit float to 16-bit half-precision float
float16_size_bytes = 16 / 8
velocity_optimized_bytes = velocity_var_count * float16_size_bytes

# 2. Optimize Color: Reduce from 32-bit float to 8-bit unsigned integer per channel
uint8_size_bytes = 8 / 8
color_optimized_bytes = color_var_count * uint8_size_bytes

# 3. Optimize Other Data: Assume it's numerical data that can be reduced from 32-bit to 16-bit floats
other_optimized_bytes = other_initial_bytes / 2

# --- Final Calculation ---
total_optimized_bytes = velocity_optimized_bytes + color_optimized_bytes + other_optimized_bytes

# --- Output the results ---
print("Calculating the optimized memory consumption per voxel:")
print(f"Initial state: {total_initial_bytes} bytes per voxel.")
print("-" * 30)

print("Optimization plan:")
print(f"1. Velocity: {int(velocity_initial_bytes)} bytes -> {int(velocity_optimized_bytes)} bytes (by using 16-bit floats).")
print(f"2. Color: {int(color_initial_bytes)} bytes -> {int(color_optimized_bytes)} bytes (by using 8-bit integers).")
print(f"3. Other Data: {int(other_initial_bytes)} bytes -> {int(other_optimized_bytes)} bytes (by halving precision).")
print("-" * 30)

# The final equation as requested
print("Final Equation (Bytes):")
print(f"{int(velocity_optimized_bytes)} (Velocity) + {int(color_optimized_bytes)} (Color) + {int(other_optimized_bytes)} (Other Data)")

# The final answer
print("\nResulting memory consumption per voxel:")
print(f"{int(total_optimized_bytes)} bytes")