def calculate_optimized_memory():
    """
    Calculates the memory consumption per voxel for an optimized smoke simulation scheme.
    """
    # 1. Define the full precision memory layout based on the problem description.
    total_bytes_full_precision = 84
    
    # Each 32-bit float is 4 bytes.
    bytes_per_float32 = 4
    
    # Velocity is twelve 32-bit floats.
    velocity_components = 12
    velocity_bytes_full = velocity_components * bytes_per_float32
    
    # Color (RGB) is three 32-bit floats.
    color_components = 3
    color_bytes_full = color_components * bytes_per_float32
    
    # Calculate the size of any other unspecified data.
    other_data_bytes_full = total_bytes_full_precision - velocity_bytes_full - color_bytes_full
    
    # 2. Apply standard optimizations.
    # Velocity and other data are optimized from 32-bit floats (4 bytes) to 16-bit halfs (2 bytes).
    bytes_per_half = 2
    velocity_bytes_optimized = velocity_components * bytes_per_half
    # We assume 'other' data is also optimizable in the same way (e.g., floats -> halfs).
    other_data_bytes_optimized = int(other_data_bytes_full / 2)
    
    # Color data is optimized from 32-bit floats (4 bytes) to 8-bit integers (1 byte).
    bytes_per_uint8 = 1
    color_bytes_optimized = color_components * bytes_per_uint8
    
    # 3. Sum the optimized components to get the final result.
    total_bytes_optimized = velocity_bytes_optimized + color_bytes_optimized + other_data_bytes_optimized
    
    print("Optimized memory consumption calculation:")
    print("The final memory per voxel is the sum of the optimized components (Velocity + Color + Other Data).")
    print("\nFinal Equation:")
    print(f"{velocity_bytes_optimized} + {color_bytes_optimized} + {other_data_bytes_optimized} = {total_bytes_optimized}")

calculate_optimized_memory()