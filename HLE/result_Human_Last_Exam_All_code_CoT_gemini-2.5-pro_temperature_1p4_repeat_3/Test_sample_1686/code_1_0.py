def solve_memory_optimization():
    """
    Calculates the optimized memory consumption per voxel in a smoke simulation.
    """
    # Step 1: Deconstruct the initial memory layout.
    initial_total_bytes = 84
    bytes_per_32bit_float = 4

    # Velocity is stored as twelve 32-bit floats.
    initial_velocity_floats = 12
    initial_velocity_bytes = initial_velocity_floats * bytes_per_32bit_float

    # Color is stored as three 32-bit floats.
    initial_color_floats = 3
    initial_color_bytes = initial_color_floats * bytes_per_32bit_float

    # The remaining bytes are for other data (e.g., density, temperature).
    initial_other_data_bytes = initial_total_bytes - initial_velocity_bytes - initial_color_bytes
    # Assuming other data is also stored as 32-bit floats.
    initial_other_data_floats = initial_other_data_bytes / bytes_per_32bit_float

    print("Initial Memory Breakdown per Voxel:")
    print(f"Total              : {initial_total_bytes} bytes")
    print(f"- Velocity         : {initial_velocity_floats} floats * {bytes_per_32bit_float} bytes/float = {initial_velocity_bytes} bytes")
    print(f"- Color            : {initial_color_floats} floats * {bytes_per_32bit_float} bytes/float = {initial_color_bytes} bytes")
    print(f"- Other Data       : {initial_other_data_bytes} bytes (inferred)")
    print("-" * 30)

    # Step 2: Define and apply the optimization strategy.
    bytes_per_16bit_half = 2  # Half-precision float
    bytes_per_8bit_char = 1   # 8-bit integer for color channels

    # Optimize velocity: Convert 32-bit floats to 16-bit half-precision floats.
    optimized_velocity_bytes = initial_velocity_floats * bytes_per_16bit_half

    # Optimize color: Convert 32-bit floats to 8-bit integers per channel.
    optimized_color_bytes = initial_color_floats * bytes_per_8bit_char

    # Optimize other data: Convert 32-bit floats to 16-bit half-precision floats.
    optimized_other_data_bytes = int(initial_other_data_floats * bytes_per_16bit_half)
    
    print("Optimized Memory Breakdown per Voxel:")
    print(f"- Velocity         : {initial_velocity_floats} floats * {bytes_per_16bit_half} bytes/half = {optimized_velocity_bytes} bytes")
    print(f"- Color            : {initial_color_floats} channels * {bytes_per_8bit_char} byte/char = {optimized_color_bytes} bytes")
    print(f"- Other Data       : {int(initial_other_data_floats)} floats * {bytes_per_16bit_half} bytes/half = {optimized_other_data_bytes} bytes")
    print("-" * 30)

    # Step 3: Calculate the final memory consumption.
    final_total_bytes = optimized_velocity_bytes + optimized_color_bytes + optimized_other_data_bytes
    
    print("Final Calculation:")
    print(f"Optimized Total = (Velocity) {optimized_velocity_bytes} + (Color) {optimized_color_bytes} + (Other Data) {optimized_other_data_bytes}")
    print(f"Resulting memory consumption per voxel: {final_total_bytes} bytes")

    return final_total_bytes

# Run the calculation and print the final answer in the required format.
final_answer = solve_memory_optimization()
print(f"<<<{final_answer}>>>")
