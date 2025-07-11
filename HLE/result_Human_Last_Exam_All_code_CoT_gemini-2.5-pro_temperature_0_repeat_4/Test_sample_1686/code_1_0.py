def calculate_optimized_memory():
    """
    Calculates the memory consumption per voxel for an optimized smoke simulation scheme.
    """
    # --- Initial Full Precision Analysis ---
    bits_in_a_byte = 8
    
    # Velocity: 12 floats * 32 bits/float
    velocity_components = 12
    initial_velocity_bits_per_component = 32
    initial_velocity_bytes = velocity_components * (initial_velocity_bits_per_component / bits_in_a_byte)
    
    # Color: 3 channels * 32 bits/channel
    color_components = 3
    initial_color_bits_per_component = 32
    initial_color_bytes = color_components * (initial_color_bits_per_component / bits_in_a_byte)
    
    # The problem states the total is 84 bytes. We calculate the unspecified "other" data.
    total_initial_bytes = 84
    other_data_bytes = total_initial_bytes - (initial_velocity_bytes + initial_color_bytes)

    # --- Optimized Scheme Calculation ---
    
    # 1. Optimize Velocity: Reduce from 32-bit float to 16-bit float (half precision)
    optimized_velocity_bits_per_component = 16
    optimized_velocity_bytes = velocity_components * (optimized_velocity_bits_per_component / bits_in_a_byte)
    
    # 2. Optimize Color: Reduce from 32-bit float to 8-bit integer (standard color depth)
    optimized_color_bits_per_component = 8
    optimized_color_bytes = color_components * (optimized_color_bits_per_component / bits_in_a_byte)
    
    # 3. The "other" data is assumed to be un-optimized as its purpose is unknown.
    
    # --- Final Calculation ---
    total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + other_data_bytes
    
    print("Optimized Memory Calculation:")
    print(f"Optimized Velocity (12 components * 16-bit): {int(optimized_velocity_bytes)} bytes")
    print(f"Optimized Color (3 components * 8-bit): {int(optimized_color_bytes)} bytes")
    print(f"Other Data (unchanged): {int(other_data_bytes)} bytes")
    print("-" * 20)
    
    # Print the final equation with each number
    print(f"Final Equation: {int(optimized_velocity_bytes)} + {int(optimized_color_bytes)} + {int(other_data_bytes)} = {int(total_optimized_bytes)}")
    print(f"The resulting memory consumption per voxel is {int(total_optimized_bytes)} bytes.")

calculate_optimized_memory()
<<<51>>>