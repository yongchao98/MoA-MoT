def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # Initial full-precision configuration
    initial_total_bytes = 84
    bytes_per_float32 = 4

    # Velocity: 12 floating-point numbers
    velocity_components = 12
    initial_velocity_bytes = velocity_components * bytes_per_float32

    # Color: 3 floating-point numbers for RGB
    color_channels = 3
    initial_color_bytes = color_channels * bytes_per_float32

    # The rest of the data (e.g., density, temperature)
    initial_other_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
    other_data_components = int(initial_other_bytes / bytes_per_float32)

    # --- Optimization Plan ---
    # Velocity and other data are reduced to 16-bit half-precision floats.
    bytes_per_float16 = 2
    # Color is reduced to 8-bit unsigned integers per channel.
    bytes_per_uint8 = 1

    # Calculate new memory size for each component
    optimized_velocity_bytes = velocity_components * bytes_per_float16
    optimized_color_bytes = color_channels * bytes_per_uint8
    optimized_other_bytes = other_data_components * bytes_per_float16

    # Calculate the final total memory per voxel
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

    # Print the breakdown of the final calculation
    print("Optimized Memory Calculation (Bytes per Voxel):")
    print(f"Velocity (12 * 16-bit) + Color (3 * 8-bit) + Other Data ({other_data_components} * 16-bit) = Total")
    print(f"{optimized_velocity_bytes} + {optimized_color_bytes} + {optimized_other_bytes} = {final_total_bytes}")

calculate_optimized_memory()
<<<39>>>