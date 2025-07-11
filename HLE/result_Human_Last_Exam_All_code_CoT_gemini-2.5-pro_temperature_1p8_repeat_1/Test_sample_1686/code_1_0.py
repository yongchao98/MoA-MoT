def calculate_optimized_memory():
    """
    Calculates the memory consumption per voxel after optimization.
    """
    # --- Step 1: Deconstruct the initial memory layout ---
    total_bytes_initial = 84
    
    # Constants for data types in the initial scheme
    BYTES_PER_FLOAT32 = 4  # 32 bits = 4 bytes
    
    # Voxel components as described
    velocity_floats_count = 12
    color_floats_count = 3
    
    # Calculate initial memory usage for specified components
    velocity_bytes_initial = velocity_floats_count * BYTES_PER_FLOAT32
    color_bytes_initial = color_floats_count * BYTES_PER_FLOAT32
    
    # The remaining bytes must be for other unspecified data (e.g., density, temperature)
    other_data_bytes_initial = total_bytes_initial - velocity_bytes_initial - color_bytes_initial
    
    # --- Step 2: Define optimization parameters ---
    
    # Common optimization is to reduce floats to half-precision (16-bit)
    BYTES_PER_FLOAT16 = 2  # 16 bits = 2 bytes
    
    # Common optimization for color is to use 8-bit integers per channel
    BYTES_PER_COLOR_UINT8 = 1 # 8 bits = 1 byte

    # --- Step 3: Calculate the optimized memory size for each component ---

    # Optimize velocity from 32-bit floats to 16-bit floats
    velocity_bytes_optimized = velocity_floats_count * BYTES_PER_FLOAT16

    # Optimize color from 32-bit floats to 8-bit integers
    color_bytes_optimized = color_floats_count * BYTES_PER_COLOR_UINT8

    # We assume the "other data" was also stored as floats and can be similarly optimized
    # by reducing precision from 32-bit to 16-bit, halving its size.
    other_data_bytes_optimized = other_data_bytes_initial // 2

    # --- Step 4: Sum the optimized components ---
    total_bytes_optimized = velocity_bytes_optimized + color_bytes_optimized + other_data_bytes_optimized

    print("Initial memory breakdown:")
    print(f"  - Velocity: {velocity_bytes_initial} bytes")
    print(f"  - Color: {color_bytes_initial} bytes")
    print(f"  - Other Data: {other_data_bytes_initial} bytes")
    print(f"  - Total: {total_bytes_initial} bytes\n")
    
    print("Optimized memory calculation:")
    print(f"  - Optimized Velocity (12 floats * 2 bytes/float): {velocity_bytes_optimized} bytes")
    print(f"  - Optimized Color (3 channels * 1 byte/channel): {color_bytes_optimized} bytes")
    print(f"  - Optimized Other Data ({other_data_bytes_initial} bytes / 2): {other_data_bytes_optimized} bytes\n")

    print("Final equation for optimized memory per voxel in bytes:")
    print(f"{velocity_bytes_optimized} + {color_bytes_optimized} + {other_data_bytes_optimized} = {total_bytes_optimized}")

calculate_optimized_memory()