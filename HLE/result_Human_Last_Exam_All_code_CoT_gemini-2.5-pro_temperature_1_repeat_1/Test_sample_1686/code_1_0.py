def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # Step 1: Deconstruct the initial memory layout.
    initial_total_bytes = 84
    float_32bit_size_bytes = 4  # 32 bits = 4 bytes
    
    # Velocity is stored as twelve 32-bit floating-point numbers.
    velocity_components = 12
    initial_velocity_bytes = velocity_components * float_32bit_size_bytes
    
    # The initial color memory is inferred from the total, as the description implies
    # velocity and color are the main components. 84 total - 48 for velocity = 36 for color.
    # This also matches a plausible, if unusual, interpretation of the problem's text
    # (3 color channels * 3 floats/channel * 4 bytes/float = 36 bytes).
    initial_color_bytes = initial_total_bytes - initial_velocity_bytes
    
    print("--- Initial Memory Layout ---")
    print(f"Initial Velocity Memory: {velocity_components} floats * {float_32bit_size_bytes} bytes/float = {initial_velocity_bytes} bytes")
    print(f"Initial Color Memory (inferred): {initial_total_bytes} total bytes - {initial_velocity_bytes} velocity bytes = {initial_color_bytes} bytes")
    print(f"Initial Total per Voxel: {initial_velocity_bytes} + {initial_color_bytes} = {initial_total_bytes} bytes")
    print("-" * 30)

    # Step 2: Apply standard optimization techniques.
    
    # For velocity, 16-bit half-precision floats are a common optimization.
    half_float_16bit_size_bytes = 2 # 16 bits = 2 bytes
    optimized_velocity_bytes = velocity_components * half_float_16bit_size_bytes
    
    # For color, 8-bit integers per channel (R, G, B) are standard for visualization.
    color_channels = 3 # R, G, B
    uint_8bit_size_bytes = 1 # 8 bits = 1 byte
    optimized_color_bytes = color_channels * uint_8bit_size_bytes
    
    print("--- Optimized Memory Layout ---")
    print(f"Optimized Velocity: Using 16-bit half-precision floats -> {velocity_components} half-floats * {half_float_16bit_size_bytes} bytes/half-float = {optimized_velocity_bytes} bytes")
    print(f"Optimized Color: Using 8-bit integers per channel -> {color_channels} channels * {uint_8bit_size_bytes} byte/channel = {optimized_color_bytes} bytes")
    print("-" * 30)

    # Step 3: Calculate the final total memory consumption.
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes
    
    print("--- Final Calculation ---")
    print(f"Resulting Memory per Voxel = Optimized Velocity Memory + Optimized Color Memory")
    # The final equation with each number explicitly stated
    print(f"Resulting Memory per Voxel = {optimized_velocity_bytes} bytes + {optimized_color_bytes} bytes = {final_total_bytes} bytes")

calculate_optimized_memory()