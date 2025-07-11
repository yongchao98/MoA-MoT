import math

def calculate_optimized_memory():
    """
    Calculates the memory consumption per voxel for an optimized smoke simulation scheme.
    """
    # --- Step 1: Analyze the Full Precision Scheme ---
    # Data sizes in bytes
    float32_size_bytes = 4

    # Components in the full precision scheme
    num_velocity_components = 12
    num_color_components = 3
    total_memory_full_bytes = 84

    # Calculate memory of known components in the full precision scheme
    memory_velocity_full_bytes = num_velocity_components * float32_size_bytes
    memory_color_full_bytes = num_color_components * float32_size_bytes

    # Calculate the memory of the remaining 'other' data (e.g., density, temperature)
    memory_other_full_bytes = total_memory_full_bytes - memory_velocity_full_bytes - memory_color_full_bytes
    
    # Determine the number of 'other' components, assuming they are also 32-bit floats
    num_other_components = math.ceil(memory_other_full_bytes / float32_size_bytes)

    # --- Step 2: Propose and Apply Optimization ---
    # Optimized data sizes in bytes based on common practices
    float16_size_bytes = 2  # Half-precision float for velocity and other scalar fields
    uint8_size_bytes = 1    # 8-bit integer for each color channel

    # --- Step 3: Calculate the Optimized Memory ---
    # Calculate the new size for each component
    memory_velocity_optimized_bytes = num_velocity_components * float16_size_bytes
    memory_color_optimized_bytes = num_color_components * uint8_size_bytes
    memory_other_optimized_bytes = num_other_components * float16_size_bytes

    # --- Step 4: Sum and Present the Result ---
    total_memory_optimized_bytes = memory_velocity_optimized_bytes + memory_color_optimized_bytes + memory_other_optimized_bytes

    print("Based on standard optimization practices (reducing precision of floats to half-precision and color to 8-bit channels):")
    print("\nOptimized Memory Calculation:")
    print(f"Velocity: {num_velocity_components} components * {float16_size_bytes} bytes/component = {memory_velocity_optimized_bytes} bytes")
    print(f"Color:    {num_color_components} channels * {uint8_size_bytes} byte/channel = {memory_color_optimized_bytes} bytes")
    print(f"Other:    {num_other_components} fields * {float16_size_bytes} bytes/field = {memory_other_optimized_bytes} bytes")
    
    print("\nFinal Equation:")
    print(f"{memory_velocity_optimized_bytes} + {memory_color_optimized_bytes} + {memory_other_optimized_bytes} = {total_memory_optimized_bytes}")
    
    print(f"\nThe resulting memory consumption per voxel would be {total_memory_optimized_bytes} bytes.")
    
    # The final answer for the system to parse
    return total_memory_optimized_bytes

# Execute the function
final_answer = calculate_optimized_memory()
# The final answer is enclosed below for parsing.
# <<<39>>>