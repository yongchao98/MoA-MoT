def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # --- Initial State as defined in the problem ---
    initial_total_bytes_per_voxel = 84

    # Initial velocity component: 12 x 32-bit floats
    num_velocity_components = 12
    bytes_per_32bit_float = 4  # 32 bits / 8 bits/byte
    initial_velocity_bytes = num_velocity_components * bytes_per_32bit_float

    # Initial color component: 3 x 32-bit floats (for R, G, B)
    num_color_channels = 3
    initial_color_bytes = num_color_channels * bytes_per_32bit_float

    # The problem states the total is 84 bytes. The explicitly described components
    # (velocity + color) do not sum to 84. This means there is other,
    # unspecified data that we will assume remains unchanged.
    other_data_bytes = initial_total_bytes_per_voxel - initial_velocity_bytes - initial_color_bytes

    # --- Optimization Assumptions ---
    # Velocity is optimized to half-precision (16-bit) floats.
    bytes_per_16bit_float = 2  # 16 bits / 8 bits/byte
    optimized_velocity_bytes = num_velocity_components * bytes_per_16bit_float

    # Color is optimized to 8-bit unsigned integers per channel.
    bytes_per_8bit_uint = 1  # 8 bits / 8 bits/byte
    optimized_color_bytes = num_color_channels * bytes_per_8bit_uint

    # --- Final Calculation ---
    # The new total is the sum of the optimized components and the unchanged data.
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + other_data_bytes

    print("Optimized Memory Calculation Per Voxel")
    print("="*40)
    print(f"Initial Velocity Memory: {initial_velocity_bytes} bytes")
    print(f"Initial Color Memory:    {initial_color_bytes} bytes")
    print(f"Other Unspecified Data:  {other_data_bytes} bytes")
    print(f"Initial Total:           {initial_total_bytes_per_voxel} bytes")
    print("-"*40)
    print("Applying Optimizations:")
    print(f"Optimized Velocity (16-bit float): {optimized_velocity_bytes} bytes")
    print(f"Optimized Color (8-bit uint):      {optimized_color_bytes} bytes")
    print(f"Unchanged Other Data:              {other_data_bytes} bytes")
    print("="*40)
    print("Final Equation for Resulting Memory Consumption:")
    print(f"{optimized_velocity_bytes} (velocity) + {optimized_color_bytes} (color) + {other_data_bytes} (other) = {final_total_bytes} bytes")

calculate_optimized_memory()
<<<51>>>