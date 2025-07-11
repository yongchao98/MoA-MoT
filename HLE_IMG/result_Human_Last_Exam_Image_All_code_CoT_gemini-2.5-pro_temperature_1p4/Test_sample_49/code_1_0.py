import sympy

def solve_circuit():
    """
    This function determines and prints the cutoff frequency for the given ladder network.
    """

    # Define the symbols for resistance and capacitance
    r, C = sympy.symbols('r C', positive=True)

    # The numbers that appear in the derived formula for the cutoff frequency
    num_sqrt_val = 3
    num_minus_val = 1
    denominator_coefficient = 2

    # Construct the symbolic expression for the cutoff frequency omega_c
    numerator = sympy.sqrt(num_sqrt_val) - num_minus_val
    denominator = denominator_coefficient * r * C
    omega_c = numerator / denominator

    # Print the final equation in a formatted string, showing each component number.
    # The result is the angular cutoff frequency omega_c.
    print("The cutoff frequency at node a0 is expressed by the following equation:")
    print(f"omega_c = (sqrt({num_sqrt_val}) - {num_minus_val}) / ({denominator_coefficient} * r * C)")
    
    # For clarity, let's also print the expression using sympy's pretty print
    print("\nSymbolic expression using sympy:")
    sympy.pprint(omega_c, use_unicode=True)

if __name__ == '__main__':
    solve_circuit()
