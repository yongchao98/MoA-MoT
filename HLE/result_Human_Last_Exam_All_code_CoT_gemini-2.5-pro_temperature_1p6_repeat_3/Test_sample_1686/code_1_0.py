def calculate_optimized_voxel_size():
    """
    Calculates the memory consumption per voxel after a standard optimization.

    This function first determines the composition of the original 84-byte voxel,
    then applies a common optimization scheme, and finally calculates the new size.
    """
    # --- Initial Full Precision State (in bytes) ---
    total_bytes_full = 84
    float32_bytes = 4  # 32 bits

    # Known components
    velocity_components = 12
    color_components = 3
    
    velocity_bytes_full = velocity_components * float32_bytes
    color_bytes_full = color_components * float32_bytes
    
    # Deduce the size and count of other simulation data
    other_bytes_full = total_bytes_full - velocity_bytes_full - color_bytes_full
    other_components_count = int(other_bytes_full / float32_bytes)

    # --- Optimized Precision State (in bytes) ---
    # Velocity and other physical quantities are reduced to 16-bit floats
    float16_bytes = 2  # 16 bits
    # Color is reduced to 8-bit unsigned integers per channel
    uint8_bytes = 1    # 8 bits

    # Calculate the new size for each part
    velocity_bytes_optimized = velocity_components * float16_bytes
    color_bytes_optimized = color_components * uint8_bytes
    other_bytes_optimized = other_components_count * float16_bytes
    
    # Calculate the total new size
    total_bytes_optimized = velocity_bytes_optimized + color_bytes_optimized + other_bytes_optimized
    
    # Print the explanation and the final equation with all numbers
    print("Based on a standard optimization scheme:")
    print("- Velocity and other physical data are reduced from 32-bit to 16-bit floats.")
    print("- Color data is reduced from 32-bit floats to 8-bit integers.")
    print("\nThe final memory consumption per voxel is calculated as follows:")
    print(f"(Velocity: {velocity_components} components * {float16_bytes} bytes) + (Color: {color_components} components * {uint8_bytes} byte) + (Other: {other_components_count} components * {float16_bytes} bytes)")
    
    # Final equation with numbers only
    print(f"({velocity_components} * {float16_bytes}) + ({color_components} * {uint8_bytes}) + ({other_components_count} * {float16_bytes}) = {total_bytes_optimized} bytes")

calculate_optimized_voxel_size()
<<<39>>>