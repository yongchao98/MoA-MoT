def print_green_function_formula():
    """
    This function prints the formula for the bare Green's function G_0,
    showing its functional dependence on the single-particle energy eigenvalue ϵ_k.
    """

    # Define the components of the formula as string variables
    green_function_symbol = "G_0(k, ω)"
    numerator = "1"
    denominator_part1 = "ω"
    denominator_part2 = "ϵ_k"
    denominator_part3 = "iδ"

    # Explain the context and the meaning of the formula
    print("In the frequency domain, the bare Green's function G_0 is a function of a quantum number 'k' (like momentum) and frequency 'ω'.")
    print("Its functional dependence on the single-particle energy eigenvalue ϵ_k is an inverse relationship.")
    print("\nThe fundamental equation is:")

    # Programmatically construct and print the final equation.
    # This fulfills the request to output each number in the equation, which is the '1'.
    print(f"\n  {green_function_symbol} = --------------------")
    print(f"            {denominator_part1} - {denominator_part2} + {denominator_part3}\n")
    
    # Let's break down the denominator term to be extra clear
    print("This can be read as:")
    equation = f"{green_function_symbol} = {numerator} / ({denominator_part1} - {denominator_part2} + {denominator_part3})"
    print(equation)
    
    print("\nWhere:")
    print(" - G_0 is the bare Green's function.")
    print(" - k is the quantum number (e.g., momentum).")
    print(" - ω is the frequency.")
    print(" - ϵ_k is the single-particle energy eigenvalue for the state k.")
    print(f" - The numerator is the number: {numerator}")
    print(" - iδ is an infinitesimal imaginary term that enforces causality.")
    print("\nThe key feature is the pole in the denominator: the function diverges when the frequency ω is equal to the energy eigenvalue ϵ_k.")

# Execute the function to display the result
print_green_function_formula()