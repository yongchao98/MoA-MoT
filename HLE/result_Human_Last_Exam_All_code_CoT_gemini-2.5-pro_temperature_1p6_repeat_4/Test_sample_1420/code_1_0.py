import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent α using the epsilon expansion.
    """
    # Define the parameters based on the problem
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension for the O(n) model
    n = 1  # Number of spin components (Ising model universality class)

    print("This script calculates the specific heat exponent alpha (α).")
    print("Method: First-order epsilon (ϵ) expansion for the O(n) model.\n")
    print(f"Given parameters:")
    print(f"  - Spatial dimension, d = {d}")
    print(f"  - Upper critical dimension, d_c = {d_c}")
    print(f"  - Number of spin components, n = {n} (assuming Ising model)\n")

    # Step 1: Calculate epsilon (ϵ)
    epsilon = d_c - d
    print(f"Step 1: Calculate the expansion parameter ϵ.")
    print(f"ϵ = d_c - d")
    print(f"ϵ = {d_c} - {d} = {epsilon}\n")

    # Step 2: Apply the formula for α
    # The formula is α = (4 - n) / (2 * (n + 8)) * ϵ
    # Using SymPy for a precise fractional representation
    n_sym = sympy.S(n)
    epsilon_sym = sympy.S(epsilon)
    
    numerator_expr = 4 - n_sym
    denominator_expr = 2 * (n_sym + 8)
    
    alpha_expr = (numerator_expr / denominator_expr) * epsilon_sym

    numerator_val = 4 - n
    denominator_val = 2 * (n + 8)
    
    print("Step 2: Use the formula α = (4 - n) / (2 * (n + 8)) * ϵ.")
    print("Substituting the values:\n")
    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = {numerator_val} / (2 * {n + 8}) * {epsilon}")
    print(f"α = {numerator_val} / {denominator_val} * {epsilon}\n")

    # Step 3: Compute the final result
    alpha_fraction = sympy.simplify(alpha_expr)
    alpha_decimal = float(alpha_fraction)

    print("Step 3: Calculate the final value.")
    print(f"The result for α as a fraction is: {alpha_fraction}")
    print(f"The result for α as a decimal is approximately: {alpha_decimal:.4f}\n")
    print("This value is the first-order approximation in the epsilon expansion.")

if __name__ == '__main__':
    calculate_alpha_exponent()