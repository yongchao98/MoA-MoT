import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent α using the first-order
    epsilon expansion for a system in d=3 dimensions.
    """
    # Step 1: Define the key parameters.
    # d is the spatial dimension.
    d = 3
    # d_c is the upper critical dimension for the n-vector model.
    d_c = 4
    # n is the number of components of the order parameter.
    # We assume the Ising model, where n=1.
    n = 1

    # Step 2: Calculate epsilon (ϵ).
    epsilon = d_c - d

    # Step 3 & 4: Use the formula and compute alpha.
    # The first-order expansion for α is: α ≈ (4 - n) / (2 * (n + 8)) * ϵ
    # We use Sympy to represent it as an exact fraction.
    numerator = (4 - n) * epsilon
    denominator = 2 * (n + 8)
    alpha_fraction = sympy.Rational(numerator, denominator)
    alpha_decimal = float(alpha_fraction)

    # Step 5: Print the explanation and the result.
    print("In renormalization group theory, the specific heat scaling exponent α can be calculated using the epsilon (ϵ) expansion.")
    print("\n--- Calculation Steps ---")
    print(f"1. The spatial dimension is d = {d}.")
    print(f"2. The upper critical dimension is d_c = {d_c}.")
    print(f"3. Epsilon is defined as ϵ = d_c - d, so ϵ = {d_c} - {d} = {epsilon}.")
    print(f"4. For the Ising model, the order parameter has n = {n} component.")
    
    print("\nThe first-order expansion formula for α is:")
    print("α ≈ (4 - n) / (2 * (n + 8)) * ϵ")

    print("\nSubstituting the values into the equation:")
    # Printing each number in the equation
    print(f"α ≈ ({4} - {n}) / ({2} * ({n} + {8})) * {epsilon}")
    
    # Showing the intermediate calculation
    print(f"α ≈ {4-n} / ({2} * {n+8}) * {epsilon}")
    print(f"α ≈ {4-n} / {2*(n+8)}")

    print(f"\nThe result is:")
    print(f"α ≈ {alpha_fraction} ≈ {alpha_decimal:.4f}")

calculate_alpha_exponent()