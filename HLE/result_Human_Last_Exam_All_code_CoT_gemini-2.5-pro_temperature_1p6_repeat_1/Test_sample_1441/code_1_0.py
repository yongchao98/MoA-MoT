def calculate_alpha_scaling(upper_critical_dim, dimensions):
    """
    Calculates and prints the specific heat critical exponent (alpha)
    based on a first-order epsilon expansion for a scalar field theory.

    The scaling relation is α ≈ (d_c - d) / 6.

    Args:
        upper_critical_dim (int): The upper critical dimension (d_c).
        dimensions (list): A list of spatial dimensions (d) to calculate for.
    """
    numerator_constant_N = 3  # (4 - N) for N=1
    denominator_constant_N = 18 # 2 * (N + 8) for N=1
    
    # The simplified fraction is 1/6
    final_denominator = 6

    print(f"In scalar field theory (Ising universality class), the specific heat exponent α")
    print(f"scales with dimensionality d (for d < {upper_critical_dim}) according to the first-order epsilon expansion.")
    print("\nThe relationship is approximately:")
    print(f"  α ≈ ({upper_critical_dim} - d) / {final_denominator}")
    print("-" * 30)

    for d in dimensions:
        if d >= upper_critical_dim:
            print(f"For d = {d}, the model is at or above the upper critical dimension. Mean-field theory applies.")
            continue
            
        # Calculate epsilon = d_c - d
        epsilon = upper_critical_dim - d
        
        # Calculate alpha ≈ epsilon / 6
        alpha = epsilon / final_denominator
        
        print(f"For spatial dimension d = {d}:")
        # Print the calculation step-by-step as requested
        print(f"  α ≈ ({upper_critical_dim} - {d}) / {final_denominator} = {epsilon} / {final_denominator} ≈ {alpha:.4f}")

if __name__ == "__main__":
    # Define the upper critical dimension for scalar phi^4 theory
    d_c = 4
    
    # Define the spatial dimensions below d_c to investigate
    dims_to_check = [1, 2, 3]
    
    # Run the calculation and display the results
    calculate_alpha_scaling(d_c, dims_to_check)
