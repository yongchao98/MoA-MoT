def calculate_alpha_scaling():
    """
    Calculates and displays the scaling relation for the specific heat
    critical exponent alpha with spatial dimensionality d.
    """
    # For a standard scalar field theory (Ising universality class), the
    # number of components of the order parameter, N, is 1.
    N = 1

    # The upper critical dimension for this model is 4.
    d_c = 4

    # The epsilon expansion gives alpha to first order in epsilon = (d_c - d).
    # The formula is: alpha = ( (4 - N) / (2 * (N + 8)) ) * epsilon

    # Let's calculate the numerator and denominator of the coefficient.
    numerator_coeff = 4 - N
    denominator_coeff = 2 * (N + 8)

    # To simplify the final expression, we find the final denominator
    # for the form: alpha = (d_c - d) / final_denominator
    final_denominator = denominator_coeff / numerator_coeff

    print("Deriving the scaling relation for the specific heat critical exponent α:")
    print(f"The calculation is for an O(N) model with N = {N}.")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print("\nThe general first-order formula in epsilon (ε = d_c - d) is:")
    print(f"α ≈ ( (4 - N) / (2 * (N + 8)) ) * ε")
    print("\nSubstituting the values N=1 and d_c=4:")
    print(f"α ≈ ( ({4} - {N}) / (2 * ({N} + {8})) ) * ({d_c} - d)")
    print(f"α ≈ ( {numerator_coeff} / {denominator_coeff} ) * ({d_c} - d)")
    print("\nThis gives the final quantitative scaling relation:")
    # The final print statement shows the equation with all numbers.
    print(f"α ≈ ({d_c} - d) / {final_denominator}")

calculate_alpha_scaling()