def calculate_optimized_memory():
    """
    Calculates the optimized memory consumption for a smoke simulation voxel.
    """
    # The original memory consumption per voxel in a full precision scheme.
    original_bytes_per_voxel = 84

    # The optimization strategy is to convert 32-bit floating-point numbers
    # to 16-bit half-precision floats. This reduces the memory size by a factor of 2.
    optimization_factor = 2

    # Calculate the new memory size after optimization.
    optimized_bytes_per_voxel = original_bytes_per_voxel // optimization_factor

    print("The original memory consumption is 84 bytes per voxel using 32-bit floats.")
    print("A standard optimization is to convert these to 16-bit half-precision floats,")
    print("which cuts the memory requirement in half while maintaining sufficient precision for simulation.")
    print("\nThe calculation for the new memory consumption is:")
    
    # Print the final equation with all the numbers
    print(f"{original_bytes_per_voxel} bytes / {optimization_factor} = {optimized_bytes_per_voxel} bytes")

    print(f"\nThe resulting memory consumption per voxel is {optimized_bytes_per_voxel} bytes.")


if __name__ == "__main__":
    calculate_optimized_memory()