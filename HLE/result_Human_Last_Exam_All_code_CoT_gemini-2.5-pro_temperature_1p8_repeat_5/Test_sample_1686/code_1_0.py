def solve_memory_optimization():
    """
    Calculates the optimized memory consumption for a smoke simulation voxel.
    """
    # Step 1: Define the initial full-precision memory layout in bytes.
    # A 32-bit float or integer is 4 bytes.
    bytes_per_32bit_float = 4
    
    # Components from the problem description
    num_velocity_floats = 12
    num_color_floats = 3
    initial_total_bytes = 84

    # Calculate memory for each component in the initial scheme
    initial_velocity_bytes = num_velocity_floats * bytes_per_32bit_float
    initial_color_bytes = num_color_floats * bytes_per_32bit_float
    
    # Infer the memory used by other unspecified fields (e.g., density, temperature)
    initial_other_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
    num_other_floats = initial_other_bytes // bytes_per_32bit_float

    # Step 2: Define the sizes for the optimized data types in bytes.
    # Half-precision float is 16 bits (2 bytes).
    # An 8-bit unsigned integer is 1 byte.
    bytes_per_16bit_float = 2
    bytes_per_8bit_integer = 1

    # Step 3: Calculate the new memory size for each optimized component.
    # Velocity and other fields are reduced to 16-bit floats.
    optimized_velocity_bytes = num_velocity_floats * bytes_per_16bit_float
    optimized_other_bytes = num_other_floats * bytes_per_16bit_float

    # Color is reduced to three 8-bit integers (one for each RGB channel).
    num_color_channels = 3
    optimized_color_bytes = num_color_channels * bytes_per_8bit_integer

    # Calculate the total optimized memory consumption
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

    # Print the final equation showing each component's contribution
    print("Based on a standard optimization scheme, the final memory consumption is calculated as follows:")
    print(f"{optimized_velocity_bytes} bytes (Velocity) + {optimized_color_bytes} bytes (Color) + {optimized_other_bytes} bytes (Other Fields) = {final_total_bytes} bytes")

solve_memory_optimization()