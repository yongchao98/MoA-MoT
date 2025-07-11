def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # Initial full-precision configuration
    initial_total_bytes = 84
    velocity_float_count = 12
    color_float_count = 3
    bytes_per_32bit_float = 4  # 32 bits = 4 bytes

    # --- Step 1 & 2: Deconstruct Initial Memory ---
    initial_velocity_bytes = velocity_float_count * bytes_per_32bit_float
    initial_color_bytes = color_float_count * bytes_per_32bit_float
    
    # The remaining bytes are for other data (e.g., density, temperature) that we assume won't be changed.
    other_data_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes

    print(f"Initial Memory Breakdown:")
    print(f"  - Velocity: {velocity_float_count} x {bytes_per_32bit_float} bytes = {initial_velocity_bytes} bytes")
    print(f"  - Color: {color_float_count} x {bytes_per_32bit_float} bytes = {initial_color_bytes} bytes")
    print(f"  - Other Data: {initial_total_bytes} - {initial_velocity_bytes} - {initial_color_bytes} = {other_data_bytes} bytes")
    print("-" * 30)

    # --- Step 3: Apply Optimizations ---
    # Optimize velocity: Use 16-bit half-precision floats (2 bytes)
    bytes_per_16bit_float = 2
    optimized_velocity_bytes = velocity_float_count * bytes_per_16bit_float

    # Optimize color: Use 8-bit unsigned integers (1 byte) per channel
    bytes_per_8bit_int = 1
    optimized_color_bytes = color_float_count * bytes_per_8bit_int

    print("Optimization Plan:")
    print(f"  - Velocity optimized to 16-bit floats: {velocity_float_count} x {bytes_per_16bit_float} bytes = {optimized_velocity_bytes} bytes")
    print(f"  - Color optimized to 8-bit integers: {color_float_count} x {bytes_per_8bit_int} byte = {optimized_color_bytes} bytes")
    print("-" * 30)

    # --- Step 4: Calculate New Total Memory ---
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + other_data_bytes

    print("Final Calculation:")
    # The final print statement shows the equation with all its components.
    print(f"Resulting Memory = {optimized_velocity_bytes} (velocity) + {optimized_color_bytes} (color) + {other_data_bytes} (other) = {final_total_bytes} bytes")

    return final_total_bytes

# Run the calculation and store the final answer.
final_answer = calculate_optimized_memory()

# The final answer is wrapped according to the required format.
print(f"\n<<<The resulting memory consumption per voxel would be {final_answer} bytes.>>>")
print(f"<<<{final_answer}>>>")