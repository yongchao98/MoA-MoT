def print_norm_bound_relation():
    """
    This function derives and prints the relationship between the maximum norm bound
    and the covolume for real quadratic fields.
    """
    
    # Define the numerator and denominator of the constant factor
    numerator = 1
    denominator = 2

    # Define the symbols for the mathematical quantities
    norm_bound_symbol = "k_{k,∞}"
    covolume_symbol = "V"

    # Explain the context of the formula
    print("For a real quadratic field defined by a squarefree natural number N, Q(sqrt(N)):")
    print(f"The upper bound for the norm in any ideal class, denoted by '{norm_bound_symbol}',")
    print(f"is related to the covolume of the ring of integers, denoted by '{covolume_symbol}',")
    print("by the following inequality derived from the Minkowski bound:")
    print("-" * 60)
    
    # Print the final equation, showing each number in the formula
    print(f"Equation: {norm_bound_symbol} <= ({numerator}/{denominator}) * {covolume_symbol}")
    
    print("-" * 60)

# Execute the function to display the answer
print_norm_bound_relation()
