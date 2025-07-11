def calculate_alpha(d, N):
    """
    Calculates the specific heat critical exponent 'alpha' using a first-order
    epsilon expansion for a given spatial dimensionality 'd' and number of
    field components 'N'.

    This approximation is valid for 'd' close to 4.
    """
    if d >= 4:
        print(f"For d = {d} (which is >= 4), mean-field theory applies.")
        print("The specific heat has a finite discontinuity, so the critical exponent alpha = 0.")
        return

    # Epsilon is the small parameter for the expansion
    epsilon = 4 - d

    # Numerator and denominator terms for the formula
    numerator = 4 - N
    denominator = 2 * (N + 8)

    # Calculate alpha using the formula
    alpha = (numerator / denominator) * epsilon

    # --- Output the calculation step-by-step as requested ---
    print(f"Calculating the critical exponent alpha for d = {d} and N = {N}:")
    print("-" * 50)
    print(f"The formula for alpha (to first order in epsilon) is:")
    print("alpha = (4 - N) / (2 * (N + 8)) * epsilon")
    print(f"where epsilon = 4 - d.")
    print("\nStep 1: Calculate epsilon")
    print(f"epsilon = 4 - {d} = {epsilon}")

    print("\nStep 2: Substitute values into the formula")
    print(f"alpha ≈ ({4} - {N}) / ({2} * ({N} + {8})) * {epsilon}")

    print("\nStep 3: Simplify the expression")
    print(f"alpha ≈ {numerator} / ({2} * {N+8}) * {epsilon}")
    print(f"alpha ≈ {numerator} / {denominator} * {epsilon}")

    print("\nStep 4: Final calculation")
    print(f"alpha ≈ {numerator / denominator} * {epsilon}")
    print(f"alpha ≈ {alpha}")
    print("-" * 50)


# --- Example Calculation ---
# You can change these values to explore different systems.
# For example:
# N=1: Ising model
# N=2: XY model
# N=3: Heisenberg model
spatial_dimensionality = 3
field_components = 1

calculate_alpha(spatial_dimensionality, field_components)

<<<alpha = (4-N)/(2*(N+8)) * (4-d)>>>