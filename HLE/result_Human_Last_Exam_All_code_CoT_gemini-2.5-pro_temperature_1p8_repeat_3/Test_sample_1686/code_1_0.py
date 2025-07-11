def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # --- Initial State (Full Precision) ---
    total_original_bytes = 84
    num_velocity_components = 12
    num_color_components = 3
    bytes_per_32_bit_float = 4

    # Calculate memory for known components and deduce the rest
    velocity_original_total_bytes = num_velocity_components * bytes_per_32_bit_float
    color_original_total_bytes = num_color_components * bytes_per_32_bit_float
    other_original_total_bytes = total_original_bytes - velocity_original_total_bytes - color_original_total_bytes
    
    # Assume the 'other' data is also composed of 32-bit floats
    num_other_components = other_original_total_bytes // bytes_per_32_bit_float

    # --- Optimized State ---
    # Velocity and Other Data are reduced to 16-bit half-precision floats (2 bytes)
    bytes_per_16_bit_float = 2
    # Color is reduced to 8-bit unsigned integers (1 byte)
    bytes_per_8_bit_uint = 1

    # Calculate the new memory consumption for each component
    velocity_optimized_total_bytes = num_velocity_components * bytes_per_16_bit_float
    color_optimized_total_bytes = num_color_components * bytes_per_8_bit_uint
    other_optimized_total_bytes = num_other_components * bytes_per_16_bit_float

    # Calculate the final total
    total_optimized_bytes = velocity_optimized_total_bytes + color_optimized_total_bytes + other_optimized_total_bytes

    print("Optimized Memory Calculation Breakdown:")
    print(f"Velocity:   {num_velocity_components} components * {bytes_per_16_bit_float} bytes/component = {velocity_optimized_total_bytes} bytes")
    print(f"Color:      {num_color_components} components * {bytes_per_8_bit_uint} byte/component  = {color_optimized_total_bytes} bytes")
    print(f"Other Data: {num_other_components} components * {bytes_per_16_bit_float} bytes/component = {other_optimized_total_bytes} bytes")
    print("\nResulting Memory Consumption Equation (in bytes):")
    print(f"{velocity_optimized_total_bytes} + {color_optimized_total_bytes} + {other_optimized_total_bytes} = {total_optimized_bytes}")

calculate_optimized_memory()
<<<39>>>