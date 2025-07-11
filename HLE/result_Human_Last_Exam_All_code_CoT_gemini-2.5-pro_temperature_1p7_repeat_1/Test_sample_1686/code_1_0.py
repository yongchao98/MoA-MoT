def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption per voxel for a smoke simulation.
    """
    # --- Initial State Analysis ---
    initial_total_bytes = 84
    
    # Velocity: 12 components as 32-bit (4-byte) floats
    velocity_components = 12
    float32_bytes = 4
    initial_velocity_bytes = velocity_components * float32_bytes
    
    # Color: 3 channels as 32-bit (4-byte) floats
    color_components = 3
    initial_color_bytes = color_components * float32_bytes
    
    # Other data is the remainder
    initial_other_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
    
    # --- Optimization Strategy & Calculation ---
    
    # Optimize Velocity: From 32-bit float to 16-bit half-precision float (2 bytes)
    float16_bytes = 2
    optimized_velocity_bytes = velocity_components * float16_bytes
    
    # Optimize Color: From 32-bit float to 8-bit unsigned integer (1 byte)
    uint8_bytes = 1
    optimized_color_bytes = color_components * uint8_bytes
    
    # Optimize Other Data: Assume a 50% reduction (e.g., all 32-bit -> 16-bit floats)
    optimized_other_bytes = initial_other_bytes / 2
    
    # --- Final Calculation ---
    total_optimized_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_bytes
    
    print("Based on a common optimization scheme:")
    print(f" - Velocity is reduced from 32-bit to 16-bit floats.")
    print(f" - Color is reduced from 32-bit floats to 8-bit integers.")
    print(f" - Other simulation data is assumed to be reduced by half.")
    print("\nThe new memory consumption is calculated as follows:\n")
    print("Final Equation:")
    # We use int() for clean output as the numbers are whole.
    print(f"{int(optimized_velocity_bytes)} (velocity) + {int(optimized_color_bytes)} (color) + {int(optimized_other_bytes)} (other) = {int(total_optimized_bytes)} bytes")
    
    print(f"\nThe resulting memory consumption per voxel is {int(total_optimized_bytes)} bytes.")

calculate_optimized_memory()