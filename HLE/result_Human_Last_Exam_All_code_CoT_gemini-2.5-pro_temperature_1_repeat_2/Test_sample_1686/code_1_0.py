def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # Initial full-precision parameters
    total_mem_initial_bytes = 84
    velocity_components = 12
    color_components = 3
    bits_per_float32 = 32
    bits_per_byte = 8

    # --- Step 1: Deconstruct the initial memory layout ---
    # Calculate memory for velocity (twelve 32-bit floats)
    mem_velocity_initial_bytes = (velocity_components * bits_per_float32) / bits_per_byte

    # Calculate memory for color (three 32-bit floats)
    mem_color_initial_bytes = (color_components * bits_per_float32) / bits_per_byte

    # The remaining memory is for other simulation data (e.g., density, temperature)
    mem_other_initial_bytes = total_mem_initial_bytes - (mem_velocity_initial_bytes + mem_color_initial_bytes)

    # --- Step 2: Define optimization parameters ---
    bits_per_float16 = 16  # Half-precision float for velocity and other data
    bits_per_uint8 = 8     # 8-bit unsigned integer for color channels

    # --- Step 3: Calculate the new, optimized memory consumption ---
    # Optimized velocity memory (twelve 16-bit floats)
    mem_velocity_optimized_bytes = (velocity_components * bits_per_float16) / bits_per_byte

    # Optimized color memory (three 8-bit integers)
    mem_color_optimized_bytes = (color_components * bits_per_uint8) / bits_per_byte
    
    # Optimized other data memory (halved by using 16-bit floats)
    mem_other_optimized_bytes = mem_other_initial_bytes / 2
    
    # --- Step 4: Calculate the final total ---
    total_mem_optimized_bytes = (
        mem_velocity_optimized_bytes +
        mem_color_optimized_bytes +
        mem_other_optimized_bytes
    )

    # --- Step 5: Print the breakdown and final result ---
    print("Optimization Calculation Breakdown:")
    print(f"Optimized Velocity: {int(mem_velocity_optimized_bytes)} bytes")
    print(f"Optimized Color: {int(mem_color_optimized_bytes)} bytes")
    print(f"Optimized Other Data: {int(mem_other_optimized_bytes)} bytes")
    print("\nFinal Equation:")
    # The prompt requires printing each number in the final equation.
    final_equation = (
        f"{int(mem_velocity_optimized_bytes)} + "
        f"{int(mem_color_optimized_bytes)} + "
        f"{int(mem_other_optimized_bytes)} = "
        f"{int(total_mem_optimized_bytes)}"
    )
    print(final_equation)
    print(f"\nResulting memory consumption per voxel: {int(total_mem_optimized_bytes)} bytes")


calculate_optimized_memory()
<<<39>>>