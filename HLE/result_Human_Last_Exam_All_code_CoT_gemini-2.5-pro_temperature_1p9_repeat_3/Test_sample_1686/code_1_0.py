def calculate_optimized_voxel_memory():
    """
    Calculates the memory consumption per voxel after a standard optimization.
    """
    # --- Component Counts ---
    # As described in the problem for the data structure to be optimized.
    num_velocity_components = 12
    num_color_components = 3

    # --- Proposed Optimized Data Types (in bytes) ---
    # Velocity is optimized from 32-bit float (4 bytes) to 16-bit half-float (2 bytes).
    # This is a common trade-off for precision vs. memory in graphics.
    bytes_per_optimized_velocity_component = 2

    # Color is optimized from 32-bit float (4 bytes) to 8-bit unsigned integer (1 byte).
    # This is standard for representing RGB color channels (0-255).
    bytes_per_optimized_color_component = 1

    # --- Calculate Memory for Each Optimized Component ---
    optimized_velocity_memory = num_velocity_components * bytes_per_optimized_velocity_component
    optimized_color_memory = num_color_components * bytes_per_optimized_color_component

    # --- Calculate Total Optimized Memory ---
    total_optimized_memory = optimized_velocity_memory + optimized_color_memory

    # --- Output the Final Equation ---
    # The final output prints each number used in the total memory calculation.
    print("Optimized Memory Calculation Per Voxel:")
    print("-------------------------------------")
    print(f"Velocity Memory: {num_velocity_components} components * {bytes_per_optimized_velocity_component} bytes/component = {optimized_velocity_memory} bytes")
    print(f"Color Memory:    {num_color_components} channels * {bytes_per_optimized_color_component} byte/channel = {optimized_color_memory} bytes")
    print("-------------------------------------")
    print(f"Total Memory = ({num_velocity_components} * {bytes_per_optimized_velocity_component}) + ({num_color_components} * {bytes_per_optimized_color_component}) = {total_optimized_memory} bytes")


calculate_optimized_voxel_memory()