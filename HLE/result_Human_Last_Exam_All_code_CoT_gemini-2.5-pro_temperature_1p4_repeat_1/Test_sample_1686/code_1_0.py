def calculate_optimized_voxel_size():
    """
    Calculates the memory consumption per voxel after optimization.
    """
    # --- Initial Full Precision State ---
    initial_total_bytes = 84
    num_velocity_components = 12
    num_color_components = 3
    bytes_per_float32 = 4 # 32 bits

    # Calculate memory for each component in the initial scheme
    initial_velocity_bytes = num_velocity_components * bytes_per_float32
    initial_color_bytes = num_color_components * bytes_per_float32
    
    # The problem implies other data exists. Let's calculate its size.
    other_data_bytes = initial_total_bytes - (initial_velocity_bytes + initial_color_bytes)

    print(f"Initial Analysis:")
    print(f"  - Velocity Memory (12 * FP32): {initial_velocity_bytes} bytes")
    print(f"  - Color Memory (3 * FP32): {initial_color_bytes} bytes")
    print(f"  - Other Unspecified Data: {initial_total_bytes} - ({initial_velocity_bytes} + {initial_color_bytes}) = {other_data_bytes} bytes")
    print("-" * 30)

    # --- Optimized State ---
    # Optimization assumptions:
    # 1. Velocity is quantized from 32-bit floats to 16-bit half-precision floats.
    # 2. Color is quantized from 32-bit floats to 8-bit unsigned integers per channel.
    # 3. Other data remains unchanged as its purpose is unknown.
    bytes_per_float16 = 2 # 16 bits
    bytes_per_uint8 = 1   # 8 bits

    # Calculate memory for each component in the optimized scheme
    optimized_velocity_bytes = num_velocity_components * bytes_per_float16
    optimized_color_bytes = num_color_components * bytes_per_uint8
    
    # The other data remains un-optimized
    final_other_data_bytes = other_data_bytes

    # Calculate the final total memory consumption
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + final_other_data_bytes

    print("Proposed Optimization:")
    print(f"  - Velocity: Quantize to 16-bit floats ({bytes_per_float16} bytes each)")
    print(f"  - Color: Quantize to 8-bit integers ({bytes_per_uint8} byte each)")
    print(f"  - Other Data: Unchanged")
    print("-" * 30)
    
    print("Resulting Memory Consumption per Voxel:")
    print("Optimized Velocity + Optimized Color + Other Data = Total")
    # The final print statement showing the full equation as requested
    print(f"{optimized_velocity_bytes} bytes + {optimized_color_bytes} bytes + {final_other_data_bytes} bytes = {final_total_bytes} bytes")

calculate_optimized_voxel_size()