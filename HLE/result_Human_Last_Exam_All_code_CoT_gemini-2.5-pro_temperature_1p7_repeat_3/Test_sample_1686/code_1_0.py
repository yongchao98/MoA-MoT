def calculate_optimized_voxel_size():
    """
    Calculates the memory consumption per voxel for an optimized smoke simulation.
    """
    # --- Step 1: Deconstruct the initial memory layout ---
    # A 32-bit float is 4 bytes.
    bytes_per_float32 = 4

    # The problem specifies velocity and color components.
    num_velocity_components = 12
    num_color_components = 3

    # The remaining bytes are assumed to be for other scalar fields.
    # Total memory = 84 bytes.
    # Specified memory = (12 * 4) + (3 * 4) = 48 + 12 = 60 bytes.
    # Unspecified memory = 84 - 60 = 24 bytes.
    # Number of other components = 24 bytes / 4 bytes/component = 6.
    num_other_components = 6

    # --- Step 2: Define optimization data types ---
    # Velocity and other fields are reduced to 16-bit half-precision floats (2 bytes).
    bytes_per_half_float = 2
    # Color is reduced to 8-bit unsigned integers (1 byte).
    bytes_per_uint8 = 1

    # --- Step 3: Calculate the new optimized memory for each part ---
    optimized_velocity_bytes = num_velocity_components * bytes_per_half_float
    optimized_color_bytes = num_color_components * bytes_per_uint8
    optimized_other_bytes = num_other_components * bytes_per_half_float

    # --- Step 4: Sum the parts and print the final equation ---
    total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes

    print("Calculating the optimized memory per voxel:")
    print(f"Optimized Velocity: {num_velocity_components} components * {bytes_per_half_float} bytes = {optimized_velocity_bytes} bytes")
    print(f"Optimized Color: {num_color_components} components * {bytes_per_uint8} byte = {optimized_color_bytes} bytes")
    print(f"Optimized Other Data: {num_other_components} components * {bytes_per_half_float} bytes = {optimized_other_bytes} bytes")
    print("-" * 40)
    print("Final equation for memory consumption per voxel:")
    # The final equation with each number explicitly shown
    print(f"{optimized_velocity_bytes} + {optimized_color_bytes} + {optimized_other_bytes} = {total_optimized_bytes}")
    print(f"\nThe resulting memory consumption per voxel is {total_optimized_bytes} bytes.")


calculate_optimized_voxel_size()