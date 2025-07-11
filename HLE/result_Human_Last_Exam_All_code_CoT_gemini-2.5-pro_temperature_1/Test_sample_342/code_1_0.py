def calculate_and_print_Dg():
    """
    Calculates and prints the first 4 terms of the sequence D_g,
    showing the details of the calculation for each term.
    """
    sequence = []
    print("This script calculates the first 4 terms of the sequence D_g.\n")

    for g in range(1, 5):
        # Calculate the product term: product_{i=1 to g} (2^(2i) - 1)
        prod_factors = []
        prod_term = 1
        for i in range(1, g + 1):
            factor = (2**(2 * i) - 1)
            prod_factors.append(factor)
            prod_term *= factor

        # Calculate the power of 2 term: 2^(g^2 + 2g)
        exponent = g**2 + 2 * g
        power_of_2 = 2**exponent
        
        # Calculate D_g
        # Use integer arithmetic throughout to maintain precision
        dg_value = power_of_2 * prod_term
        sequence.append(dg_value)
        
        # Build the strings for printing the equation
        prod_factors_str = " * ".join(map(str, prod_factors))
        
        print(f"For g = {g}:")
        print(f"D_{g} = 2^({g}^2 + 2*{g}) * product_{{i=1..{g}}}(2^(2i) - 1)")
        print(f"D_{g} = 2^({exponent}) * ({prod_factors_str})")
        print(f"D_{g} = {power_of_2} * {prod_term}")
        print(f"D_{g} = {dg_value}\n")

    return sequence

# Run the calculation and store the final sequence
final_sequence = calculate_and_print_Dg()

# Format the final answer as requested
final_answer_str = ", ".join(map(str, final_sequence))
print(f"The sequence of the first 4 terms of D_g is: {final_answer_str}")

print(f"\n<<<{final_answer_str}>>>")