def calculate_optimized_voxel_size():
    """
    Calculates the memory consumption per voxel for an optimized smoke simulation.
    """
    # 1. Deconstruct the initial full-precision scheme.
    initial_total_bytes = 84

    # Velocity: 12 variables * 4 bytes/variable (32-bit float)
    num_velocity_vars = 12
    bytes_per_velocity_var_initial = 4
    initial_velocity_bytes = num_velocity_vars * bytes_per_velocity_var_initial

    # Color: 3 channels * 4 bytes/channel (32-bit float)
    num_color_vars = 3
    bytes_per_color_var_initial = 4
    initial_color_bytes = num_color_vars * bytes_per_color_var_initial

    # The prompt specifies data that sums to less than the total.
    # We assume the remaining data is for other properties (e.g., density, temperature)
    # and will not be optimized.
    initial_known_components_bytes = initial_velocity_bytes + initial_color_bytes
    other_data_bytes = initial_total_bytes - initial_known_components_bytes

    # 2. Define the optimized scheme.
    # Velocity is optimized to 16-bit half-precision floats (2 bytes).
    bytes_per_velocity_var_optimized = 2
    # Color is optimized to 8-bit unsigned integers (1 byte per channel).
    bytes_per_color_var_optimized = 1

    # 3. Calculate the new memory size for the optimized components.
    optimized_velocity_bytes = num_velocity_vars * bytes_per_velocity_var_optimized
    optimized_color_bytes = num_color_vars * bytes_per_color_var_optimized

    # 4. Calculate the final total memory consumption per voxel.
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + other_data_bytes
    
    # 5. Print the breakdown and the final result.
    print("Calculating the new memory consumption per voxel...")
    print("Equation: (Optimized Velocity Size) + (Optimized Color Size) + (Other Data Size)")
    print(f"Final calculation: {optimized_velocity_bytes} bytes + {optimized_color_bytes} bytes + {other_data_bytes} bytes = {final_total_bytes} bytes")

calculate_optimized_voxel_size()
<<<51>>>