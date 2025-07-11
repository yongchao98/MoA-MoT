def solve_memory_consumption():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # --- Step 1 & 2: Analyze the Initial Full Precision Layout ---

    # A 32-bit float is 4 bytes (32 / 8)
    bytes_per_32_bit_float = 4
    
    # Velocity is twelve 32-bit floats
    initial_velocity_components = 12
    initial_velocity_bytes = initial_velocity_components * bytes_per_32_bit_float
    
    # Color is three 32-bit floats
    initial_color_components = 3
    initial_color_bytes = initial_color_components * bytes_per_32_bit_float
    
    # The total memory is given as 84 bytes
    initial_total_bytes = 84
    
    # The remaining bytes are assumed to be for other data not subject to this optimization
    unspecified_data_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
    
    # --- Step 3: Define and Apply Optimizations ---
    
    # Optimize Velocity: Reduce from 32-bit (4 bytes) to 16-bit half-precision float (2 bytes)
    bytes_per_16_bit_float = 2
    optimized_velocity_bytes = initial_velocity_components * bytes_per_16_bit_float
    
    # Optimize Color: Reduce from 32-bit float (4 bytes) to 8-bit integer (1 byte) per channel
    bytes_per_8_bit_channel = 1
    optimized_color_bytes = initial_color_components * bytes_per_8_bit_channel
    
    # --- Step 4: Calculate the Final Memory Size ---
    
    # The final total is the sum of the optimized components and the unchanged data
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + unspecified_data_bytes

    # --- Final Output ---
    
    print("This script calculates the new memory per voxel after optimization.")
    print(f"The optimized scheme uses {bytes_per_16_bit_float} bytes for each of the {initial_velocity_components} velocity numbers, totaling {optimized_velocity_bytes} bytes.")
    print(f"Color is optimized to {bytes_per_8_bit_channel} byte per channel for {initial_color_components} channels, totaling {optimized_color_bytes} bytes.")
    print(f"The unspecified simulation data remains at {unspecified_data_bytes} bytes.")
    print("\nThe final memory consumption is the sum of these components:")
    
    # Print the final equation with each number explained
    print(f"{optimized_velocity_bytes} (velocity) + {optimized_color_bytes} (color) + {unspecified_data_bytes} (other data) = {final_total_bytes} bytes")

solve_memory_consumption()