def print_inner_product_formula():
    """
    Prints the formula for the inner product (ϕ, Dϕ) in finite-temperature field theory
    and explains its components.
    """

    # Define the components of the formula as strings for display
    temperature_factor = "T"
    summation = "Σ_{n=-∞ to ∞}"
    integration = "∫ d³k / (2π)³"
    propagator_inverse = "(ω_n² + |k|² + m²)"
    field_fourier_mode = "|ϕ̃(ω_n, k)|²"
    matsubara_frequency_def = "ω_n = 2πnT"

    # Assemble the final formula string
    formula = (
        f"{temperature_factor} {summation} {integration} "
        f"{propagator_inverse} * {field_fourier_mode}"
    )

    print("In finite-temperature field theory for a neutral scalar field, the inner product (ϕ, Dϕ)")
    print("is twice the free-field action, S_free[ϕ]. In the momentum-space representation,")
    print("this is given by:")
    print("\n(ϕ, Dϕ) = " + formula + "\n")

    print("Here is a breakdown of the components:")
    print("-" * 30)
    print(f"T: The temperature of the system.")
    print(f"n: An integer, indexing the Matsubara modes.")
    print(f"k: The 3-dimensional momentum vector.")
    print(f"m: The mass of the scalar field.")
    print(f"ϕ̃(ω_n, k): The Fourier-space component of the field ϕ for a given Matsubara")
    print(f"             frequency ω_n and momentum k.")
    
    print("\nThe mathematical operators are:")
    print(f"Σ_{{...}}: A summation over all integer values of 'n', from -∞ to ∞.")
    print(f"∫ d³k / (2π)³: An integral over all possible 3-momenta 'k'.")
    
    print("\nThe specific numbered terms are:")
    # Per the prompt, explicitly mention the numbers in the formula
    print("1. The '3' in d³k signifies the integral is over 3 spatial dimensions.")
    print("2. The '(2π)³' is the standard phase-space normalization factor.")
    print("   - This term contains the number '2'.")
    print(f"3. The term '{matsubara_frequency_def}' defines the Matsubara frequency.")
    print("   - This definition contains the number '2'.")
    print("-" * 30)

if __name__ == '__main__':
    print_inner_product_formula()