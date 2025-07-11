def solve_quantization_scheme():
    """
    Calculates the optimized memory consumption for a smoke simulation voxel.
    """
    # --- Initial Full Precision State ---
    initial_total_bytes = 84
    
    # Velocity: 12 floats * 4 bytes/float (32-bit)
    velocity_components = 12
    bytes_per_float32 = 4
    initial_velocity_bytes = velocity_components * bytes_per_float32
    
    # Color: 3 floats * 4 bytes/float (32-bit)
    color_components = 3
    initial_color_bytes = color_components * bytes_per_float32
    
    # Unspecified "other" data
    initial_other_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
    
    # --- Optimization Strategy ---
    # Velocity and other float data are reduced to half-precision (16-bit) floats.
    bytes_per_float16 = 2
    # Color data is reduced to 8-bit unsigned integers per channel.
    bytes_per_uint8 = 1
    
    # --- Optimized Memory Calculation ---
    
    # Optimized Velocity: 12 components * 2 bytes/half-float
    optimized_velocity_bytes = velocity_components * bytes_per_float16
    
    # Optimized Color: 3 components * 1 byte/uint8
    optimized_color_bytes = color_components * bytes_per_uint8
    
    # Optimized Other Data: Assuming the "other" data was also float-based,
    # we can reduce it to half-precision as well.
    # Number of 'other' components = 24 bytes / 4 bytes/float = 6 components
    other_data_components = initial_other_bytes // bytes_per_float32
    optimized_other_bytes = other_data_components * bytes_per_float16
    
    # Total optimized memory
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes
    
    print("Based on standard optimization techniques (using 16-bit floats for physical quantities and 8-bit integers for color channels):")
    print("\nResulting Memory Consumption per Voxel:")
    print(f"{optimized_velocity_bytes} bytes (Velocity) + {optimized_color_bytes} bytes (Color) + {optimized_other_bytes} bytes (Other Data) = {final_total_bytes} bytes")

solve_quantization_scheme()