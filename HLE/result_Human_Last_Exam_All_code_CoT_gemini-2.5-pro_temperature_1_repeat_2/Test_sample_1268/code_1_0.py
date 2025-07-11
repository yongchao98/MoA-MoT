import math

def display_upper_bound_formula():
    """
    This function explains and prints the upper bound for the maximum norm
    in relation to the covolume for lattices associated with squarefree integers.
    """

    # Define the symbols used in the formula.
    # k_{k,∞} represents the maximum norm of the shortest non-zero lattice vector.
    # V represents the covolume of the lattice.
    maximum_norm_symbol = "k_{k,∞}"
    covolume_symbol = "V"

    # The relationship is k_{k,∞} <= sqrt(V), which can be written with exponents.
    # The numbers in the equation are the 1 and 2 in the exponent 1/2.
    exponent_numerator = 1
    exponent_denominator = 2

    # Print the full formula.
    print("The upper bound for the maximum norm in relation to the covolume is given by the formula:")
    # Using Unicode for better mathematical representation
    # k_{k,∞} ≤ V^(1/2)
    print(f"{maximum_norm_symbol.replace('∞', '\u221e')} \u2264 {covolume_symbol}^({exponent_numerator}/{exponent_denominator})")

    print("\n--- Formula Components ---")
    print(f"Left-hand side (Maximum Norm): {maximum_norm_symbol.replace('∞', '\u221e')}")
    print("Inequality operator: \u2264 (less than or equal to)")
    print(f"Right-hand side base (Covolume): {covolume_symbol}")
    print(f"Exponent Numerator: {exponent_numerator}")
    print(f"Exponent Denominator: {exponent_denominator}")

# Execute the function to display the result.
display_upper_bound_formula()