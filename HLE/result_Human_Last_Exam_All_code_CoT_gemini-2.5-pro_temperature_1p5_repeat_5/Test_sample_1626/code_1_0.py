import math

def demonstrate_tanh_for_grid_spacing():
    """
    This function demonstrates how the tanh function can be used to map
    computational grid points to physical grid points, creating clustering.
    It prints the steps of a sample calculation as requested.
    """
    
    # The tanh mapping function is often of the form:
    # x_physical = L * (1 + tanh(alpha * x_computational) / tanh(alpha))
    # For simplicity, we'll demonstrate the core of the mapping:
    # y = tanh(alpha * x) where x is in [-1, 1]
    
    print("A common technique for grid spacing uses the tanh function.")
    print("Here is an example calculation for a single grid point.")
    
    # Parameters for the example equation: y = tanh(alpha * x)
    alpha = 2.5  # Clustering factor. Higher alpha means more clustering.
    x_comp = 0.5  # A sample point in the "computational" uniform grid.
    
    # Step-by-step calculation
    product = alpha * x_comp
    result = math.tanh(product)
    
    # Outputting the equation with each number, as requested.
    print("\nFinal Equation Calculation:")
    print(f"y = tanh(alpha * x)")
    print(f"Using values: alpha = {alpha}, x = {x_comp}")
    print(f"y = tanh({alpha} * {x_comp})")
    print(f"y = tanh({product})")
    print(f"y = {result}")
    print("\nThis shows how a point at 0.5 in a uniform grid can be mapped to ~0.848,")
    print("demonstrating the non-linear spacing effect of the tanh function.")

demonstrate_tanh_for_grid_spacing()