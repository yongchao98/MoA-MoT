def display_bare_greens_function_formula():
    """
    This function programmatically constructs and prints the formula for the
    bare Green's function, G_0, to demonstrate its functional dependence on
    the single-particle energy eigenvalue, epsilon_k.
    """

    # --- Define the components of the formula ---
    
    # The left side of the equation
    green_function_symbol = "G_0(k, omega)"
    
    # The numerator of the fraction. This is the only number in the equation.
    numerator = 1
    
    # The symbolic variables in the denominator
    frequency_symbol = "omega"
    energy_eigenvalue_symbol = "epsilon_k"
    infinitesimal_term_symbol = "i*delta"

    # --- Print the explanation and constructed formula ---
    
    print("In the Feynman path integral formalism, the bare Green's function G_0(k, omega)")
    print("for a non-interacting system has an inverse dependence on the single-particle")
    print(f"energy eigenvalue, {energy_eigenvalue_symbol}.")
    print("\nIn frequency-momentum space, the equation is:")
    
    # Programmatically construct and print the final equation
    # The format f-string helps build the mathematical expression from its parts.
    final_equation = f"{green_function_symbol} = {numerator} / ({frequency_symbol} - {energy_eigenvalue_symbol} + {infinitesimal_term_symbol})"
    print(f"\n{final_equation}\n")
    
    print("Here is a breakdown of the final equation's components:")
    print(f"- The number in the numerator is: {numerator}")
    print("- The term 'omega' is the frequency of the particle.")
    print(f"- The term '{energy_eigenvalue_symbol}' is the single-particle energy eigenvalue for state 'k'.")
    print(f"- The term '{infinitesimal_term_symbol}' is an infinitesimal imaginary component required for the correct time-ordering.")

# Execute the function to display the result.
display_bare_greens_function_formula()
