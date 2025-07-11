def solve_memory_optimization():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # --- Step 1: Analyze the initial unoptimized state ---

    # Given initial total memory per voxel in bytes
    initial_total_memory = 84

    # Details of the unoptimized data structure
    velocity_float_count = 12
    color_float_count = 3
    bytes_per_32bit_float = 4

    # Calculate memory for known components in the initial scheme
    initial_velocity_memory = velocity_float_count * bytes_per_32bit_float
    initial_color_memory = color_float_count * bytes_per_32bit_float

    # Infer the size of the remaining data
    initial_other_data_memory = initial_total_memory - initial_velocity_memory - initial_color_memory

    # --- Step 2: Define and apply optimization strategy ---

    # Optimization for floating-point data (Velocity and Other):
    # Convert from 32-bit floats (4 bytes) to 16-bit half-precision floats (2 bytes).
    bytes_per_16bit_float = 2
    optimized_velocity_memory = velocity_float_count * bytes_per_16bit_float
    # Assume 'other data' is also composed of floating-point values that can be halved.
    optimized_other_data_memory = initial_other_data_memory // 2

    # Optimization for color data:
    # Convert from 32-bit floats to 8-bit unsigned integers (1 byte) per channel.
    color_channel_count = 3
    bytes_per_8bit_channel = 1
    optimized_color_memory = color_channel_count * bytes_per_8bit_channel

    # --- Step 3: Calculate the final total optimized memory ---

    final_total_memory = optimized_velocity_memory + optimized_color_memory + optimized_other_data_memory

    print("Calculating optimized memory per voxel:")
    print("Based on standard optimizations (32-bit -> 16-bit floats, 32-bit float colors -> 8-bit integer colors):")
    print("\nFinal Memory Calculation:")
    print(f"{optimized_velocity_memory} bytes (velocity) + {optimized_color_memory} bytes (color) + {optimized_other_data_memory} bytes (other data) = {final_total_memory} bytes")

solve_memory_optimization()