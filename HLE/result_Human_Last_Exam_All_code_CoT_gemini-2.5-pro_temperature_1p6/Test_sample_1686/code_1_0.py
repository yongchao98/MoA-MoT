import math

# Step 1: Define the memory allocation for the full precision scheme.
initial_total_per_voxel = 84  # bytes
num_velocity_components = 12
num_color_channels = 3
bytes_per_float32 = 4  # A 32-bit float is 4 bytes.

# Calculate the initial memory used by the specified components.
initial_velocity_memory = num_velocity_components * bytes_per_float32
initial_color_memory = num_color_channels * bytes_per_float32

# Determine the size of the unmentioned data.
other_data_memory = initial_total_per_voxel - (initial_velocity_memory + initial_color_memory)

print("--- Initial Memory Breakdown ---")
print(f"Initial Velocity (12 * 32-bit float): {initial_velocity_memory} bytes")
print(f"Initial Color (3 * 32-bit float): {initial_color_memory} bytes")
print(f"Other Data (by deduction): {other_data_memory} bytes")
print("-" * 34)

# Step 2: Define the memory allocation for the optimized scheme.
# For velocity, we switch from 32-bit floats to 16-bit half-precision floats.
bytes_per_float16 = 2
# For color, we switch from 32-bit floats to 8-bit unsigned integers.
bytes_per_uint8 = 1

# Calculate the optimized memory for the components.
optimized_velocity_memory = num_velocity_components * bytes_per_float16
optimized_color_memory = num_color_channels * bytes_per_uint8

print("\n--- Optimized Memory Calculation ---")
print("Strategy: Reduce velocity to 16-bit floats and color to 8-bit integers.")
print(f"Optimized Velocity Memory (12 * 16-bit float): {optimized_velocity_memory} bytes")
print(f"Optimized Color Memory (3 * 8-bit int): {optimized_color_memory} bytes")
print(f"Other Data Memory (unchanged): {other_data_memory} bytes")
print("-" * 34)


# Step 3: Calculate the final total memory per voxel.
final_total_memory = optimized_velocity_memory + optimized_color_memory + other_data_memory

# Print the final equation as requested.
print("\nFinal Equation for Optimized Memory per Voxel:")
print(f"{optimized_velocity_memory} (Velocity) + {optimized_color_memory} (Color) + {other_data_memory} (Other Data) = {final_total_memory} bytes")

print(f"\nThe resulting memory consumption per voxel is {final_total_memory} bytes.")