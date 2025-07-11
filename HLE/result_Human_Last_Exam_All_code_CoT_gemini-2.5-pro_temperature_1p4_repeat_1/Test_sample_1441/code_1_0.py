def calculate_alpha_scaling(d):
    """
    Calculates the specific heat critical exponent alpha (α)
    for a given spatial dimensionality d < 4 using first-order epsilon expansion.
    """
    if d >= 4:
        print("This first-order epsilon expansion is valid for d < 4.")
        return

    # For a single-component scalar field, the component number N is 1.
    N = 1

    # Epsilon (ε) is the deviation from the upper critical dimension of 4.
    epsilon = 4 - d

    # The formula for α to first order in ε is: (4 - N) / (2 * (N + 8)) * ε
    numerator_N = 4 - N
    denominator_N_term = N + 8
    denominator_N = 2 * denominator_N_term

    # Calculate alpha
    alpha = (numerator_N / denominator_N) * epsilon

    # Print the step-by-step calculation
    print(f"Calculating the specific heat critical exponent α for spatial dimension d = {d}")
    print("-" * 60)
    print("The scaling relation from first-order epsilon expansion is:")
    print("α = (4 - N) / (2 * (N + 8)) * (4 - d)")
    print("\nFor a single-component scalar field, we set N = 1:")
    print(f"α = (4 - {N}) / (2 * ({N} + 8)) * (4 - {d})")
    print("\nEvaluating the terms:")
    print(f"The term (4 - d) becomes: 4 - {d} = {epsilon}")
    print(f"The term (4 - N) becomes: 4 - {N} = {numerator_N}")
    print(f"The term (N + 8) becomes: {N} + 8 = {denominator_N_term}")
    print(f"The term 2 * (N + 8) becomes: 2 * {denominator_N_term} = {denominator_N}")
    print("\nSubstituting the evaluated numbers back into the equation:")
    print(f"α = {numerator_N} / {denominator_N} * {epsilon}")
    print(f"α = {numerator_N / denominator_N} * {epsilon}")
    print(f"α = {alpha}")

# Run the calculation for d=3
calculate_alpha_scaling(3)