import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.
    """
    # Define symbolic variables
    n_sym, epsilon_sym = sympy.symbols('n epsilon')

    # Parameters for the calculation
    d = 3       # Spatial dimension
    d_c = 4     # Upper critical dimension
    n = 1       # Number of order parameter components (Ising model)

    # Calculate epsilon
    epsilon = d_c - d

    # Formula for alpha to first order in epsilon
    alpha_formula_sym = (4 - n_sym) / (2 * (n_sym + 8)) * epsilon_sym
    
    # Substitute values into the formula
    alpha_value = alpha_formula_sym.subs({n_sym: n, epsilon_sym: epsilon})
    
    # --- Output Section ---
    print("This script calculates the specific heat scaling exponent α using the epsilon (ε) expansion.")
    print("-" * 70)
    print(f"The given spatial dimension is d = {d}.")
    print(f"The upper critical dimension for this universality class is d_c = {d_c}.")
    print(f"Therefore, the expansion parameter is ε = d_c - d = {d_c} - {d} = {epsilon}.")
    print("\nFor a model with 'n' order parameter components, the general formula for α to first order is:")
    print(f"α ≈ {sympy.pretty(alpha_formula_sym, use_unicode=False)}")
    
    print(f"\nWe will use the Ising model universality class, where n = {n}.")
    print("Substituting the values of n and ε into the formula:")
    
    # Show step-by-step calculation
    step1 = f"α ≈ (4 - {n}) / (2 * ({n} + 8)) * {epsilon}"
    step2 = f"α ≈ {4 - n} / (2 * {n + 8}) * {epsilon}"
    step3 = f"α ≈ {4 - n} / {2 * (n + 8)} * {epsilon}"
    step4 = f"α ≈ {sympy.nsimplify(alpha_value)}"
    
    print(step1)
    print(step2)
    print(step3)
    print(step4)

    # Print the final result in decimal form
    print(f"\nAs a decimal, the value is approximately: {float(alpha_value):.4f}")

# Execute the function
calculate_alpha_exponent()