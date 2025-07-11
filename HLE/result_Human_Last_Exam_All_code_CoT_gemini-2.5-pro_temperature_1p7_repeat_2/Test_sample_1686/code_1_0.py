def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption for a smoke simulation voxel
    based on a standard quantization scheme.
    """
    # --- Step 1: Deconstruct the initial memory layout ---
    initial_total_bytes_per_voxel = 84
    velocity_num_fields = 12  # Given as twelve 32-bit floats
    color_num_channels = 3    # Given as three 32-bit floats for RGB
    bytes_per_float32 = 4

    # Calculate memory for the known components in the full-precision scheme
    initial_velocity_bytes = velocity_num_fields * bytes_per_float32
    initial_color_bytes = color_num_channels * bytes_per_float32

    # Deduce the size of other unspecified data fields
    initial_other_bytes = initial_total_bytes_per_voxel - initial_velocity_bytes - initial_color_bytes

    print("--- Initial Memory Layout (per Voxel) ---")
    print(f"Total given: {initial_total_bytes_per_voxel} bytes")
    print(f"Velocity: {initial_velocity_bytes} bytes ({velocity_num_fields} x 32-bit floats)")
    print(f"Color: {initial_color_bytes} bytes ({color_num_channels} x 32-bit floats)")
    print(f"Other Data: {initial_other_bytes} bytes (Deduced)")
    print("\n")


    # --- Step 2: Define and apply optimizations ---
    print("--- Applying Optimizations ---")
    print("Strategy:")
    print("- Velocity/Other Fields: Converted from 32-bit float to 16-bit half-precision float.")
    print("- Color: Converted from 32-bit float to 8-bit unsigned integer per channel.")
    print("\n")
    
    bytes_per_float16 = 2  # Half precision
    bytes_per_uint8 = 1    # 8-bit integer

    # Calculate the memory for each component in the optimized scheme
    optimized_velocity_bytes = velocity_num_fields * bytes_per_float16
    optimized_color_bytes = color_num_channels * bytes_per_uint8
    
    # Other physical data is also reduced to half-precision
    optimized_other_bytes = initial_other_bytes // 2


    # --- Step 3: Calculate the final total and present the equation ---
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

    print("--- Final Optimized Memory Calculation (per Voxel) ---")
    print("The final memory consumption is the sum of the optimized components.")
    print("Equation (Velocity + Color + Other):")
    
    # Output the final equation with each number
    print(f"{optimized_velocity_bytes} + {optimized_color_bytes} + {optimized_other_bytes} = {final_total_bytes}")


if __name__ == "__main__":
    calculate_optimized_memory()