import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the one-loop epsilon expansion.
    """
    # Step 1 & 2: Define parameters
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension
    epsilon = d_c - d  # Epsilon parameter

    # Step 3 & 4: Use the formula for the Ising model (n=1)
    n = 1  # Number of spin components for the Ising model

    # Use sympy for a fractional representation
    n_sym = sympy.Integer(n)
    epsilon_sym = sympy.Integer(epsilon)

    # One-loop formula for alpha
    alpha_formula_num = 4 - n_sym
    alpha_formula_den = 2 * (n_sym + 8)
    alpha = (alpha_formula_num / alpha_formula_den) * epsilon_sym

    # Print the explanation and calculation steps
    print("Calculating the specific heat exponent α using the epsilon expansion for d=3.")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The spatial dimension is d = {d}.")
    print(f"The expansion parameter is ϵ = d_c - d = {d_c} - {d} = {epsilon}.")
    print("\nThe one-loop formula for α is: α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("Assuming the Ising model, we set n = 1.")
    print("\nSubstituting the values into the formula:")
    print(f"α = (4 - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {alpha_formula_num} / (2 * {n+8}) * {epsilon}")
    print(f"α = {alpha_formula_num} / {alpha_formula_den} * {epsilon}")
    print(f"α = {alpha} (or approximately {float(alpha):.4f})")

calculate_alpha_exponent()