import math

# --- Step 1: Deconstruct the initial memory layout ---

# Initial parameters from the problem description
total_memory_initial_bytes = 84
num_velocity_fields = 12
num_color_fields = 3

# A standard 32-bit float uses 4 bytes
bytes_per_32bit_float = 4

# Calculate memory for each component in the initial scheme
velocity_memory_initial_bytes = num_velocity_fields * bytes_per_32bit_float
color_memory_initial_bytes = num_color_fields * bytes_per_32bit_float
other_data_memory_initial_bytes = total_memory_initial_bytes - velocity_memory_initial_bytes - color_memory_initial_bytes

print("--- Initial Memory Breakdown ---")
print(f"Initial total memory: {total_memory_initial_bytes} bytes")
print(f"Velocity ({num_velocity_fields} x 32-bit floats): {velocity_memory_initial_bytes} bytes")
print(f"Color ({num_color_fields} x 32-bit floats): {color_memory_initial_bytes} bytes")
print(f"Inferred other data: {other_data_memory_initial_bytes} bytes\n")

# --- Step 2: Apply optimizations ---

print("--- Optimization Assumptions ---")
print("To reduce memory while maintaining sufficient precision:")
print("- Physical data (velocity, etc.) is optimized from 32-bit to 16-bit floats.")
print("- Color data (RGB) is optimized to one 8-bit integer per channel.\n")

# A 16-bit float uses 2 bytes. An 8-bit integer uses 1 byte.
# Optimization for physical data (velocity, other data) is a 50% reduction.
velocity_memory_optimized_bytes = velocity_memory_initial_bytes / 2
other_data_memory_optimized_bytes = other_data_memory_initial_bytes / 2

# Optimization for color: 3 channels (RGB) * 1 byte/channel
num_color_channels = 3
bytes_per_8bit_channel = 1
color_memory_optimized_bytes = num_color_channels * bytes_per_8bit_channel


# --- Step 3: Calculate the new total memory ---
final_total_memory_bytes = velocity_memory_optimized_bytes + color_memory_optimized_bytes + other_data_memory_optimized_bytes

print("--- Final Calculation ---")
print("The resulting memory consumption per voxel is calculated by summing the optimized components:")
# The final line prints the equation with each number as requested
print(f"{int(velocity_memory_optimized_bytes)} (Velocity) + {int(color_memory_optimized_bytes)} (Color) + {int(other_data_memory_optimized_bytes)} (Other) = {int(final_total_memory_bytes)}")
