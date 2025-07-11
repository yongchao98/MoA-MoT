def calculate_optimized_voxel_size():
    """
    Calculates the memory consumption per voxel for an optimized smoke simulation scheme.
    """
    # --- Step 1: Deconstruct the initial memory layout ---
    total_initial_bytes = 84
    float32_bytes = 4  # 32-bit float

    # Initial components
    num_velocity_floats = 12
    num_color_floats = 3

    initial_velocity_bytes = num_velocity_floats * float32_bytes
    initial_color_bytes = num_color_floats * float32_bytes

    # Calculate the size of other unspecified data (e.g., density, temperature)
    initial_other_bytes = total_initial_bytes - initial_velocity_bytes - initial_color_bytes
    # Assuming the other data is also stored as 32-bit floats
    num_other_floats = int(initial_other_bytes / float32_bytes)

    # --- Step 2: Define optimized data type sizes ---
    float16_bytes = 2  # 16-bit half-precision float
    uint8_bytes = 1    # 8-bit unsigned integer for color channels

    # --- Step 3: Calculate the optimized memory for each component ---
    optimized_velocity_bytes = num_velocity_floats * float16_bytes
    optimized_color_bytes = num_color_floats * uint8_bytes
    optimized_other_bytes = num_other_floats * float16_bytes

    # --- Step 4: Sum the components for the final result ---
    total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

    # --- Final Output ---
    # The final print statement shows the full equation with the numbers used.
    print(f"Optimized Voxel Size Calculation:")
    print(f"({num_velocity_floats} velocity components * {float16_bytes} bytes/component) + "
          f"({num_color_floats} color channels * {uint8_bytes} byte/channel) + "
          f"({num_other_floats} other components * {float16_bytes} bytes/component)")
    print(f"= {optimized_velocity_bytes} bytes (velocity) + "
          f"{optimized_color_bytes} bytes (color) + "
          f"{optimized_other_bytes} bytes (other)")
    print(f"= {total_optimized_bytes} bytes per voxel")

calculate_optimized_voxel_size()