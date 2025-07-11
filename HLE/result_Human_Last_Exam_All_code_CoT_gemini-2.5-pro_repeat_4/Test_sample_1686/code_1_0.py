def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # Step 1: Deconstruct the initial memory layout in bytes.
    # A 32-bit float is 4 bytes.
    velocity_original_bytes = 12 * 4
    color_original_bytes = 3 * 4
    total_original_bytes = 84
    
    # The remaining bytes are for other data (e.g., density, temperature).
    other_data_original_bytes = total_original_bytes - velocity_original_bytes - color_original_bytes
    
    # Step 2: Apply standard optimization techniques.
    
    # Optimize velocity from 32-bit floats to 16-bit (2-byte) half-precision floats.
    velocity_optimized_bytes = 12 * 2
    
    # Optimize color from 32-bit floats to 8-bit (1-byte) unsigned integers per channel.
    color_optimized_bytes = 3 * 1
    
    # Assume other data are floats and optimize them to 16-bit half-precision as well.
    other_data_optimized_bytes = other_data_original_bytes / 2
    
    # Step 3: Calculate the total optimized memory consumption.
    total_optimized_bytes = velocity_optimized_bytes + color_optimized_bytes + other_data_optimized_bytes
    
    # Print the final equation with each component clearly labeled.
    print("Optimized memory consumption calculation:")
    print(f"{velocity_optimized_bytes} bytes (velocity) + {color_optimized_bytes} bytes (color) + {int(other_data_optimized_bytes)} bytes (other data) = {int(total_optimized_bytes)} bytes per voxel")

calculate_optimized_memory()