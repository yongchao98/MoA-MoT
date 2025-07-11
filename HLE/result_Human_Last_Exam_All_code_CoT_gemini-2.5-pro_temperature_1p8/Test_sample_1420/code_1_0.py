def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha (α) using the epsilon expansion.
    """
    # --- Parameters ---
    # Spatial dimension
    d = 3
    # Upper critical dimension (for the n-vector model)
    d_c = 4
    # Number of components of the order parameter. We assume the Ising model.
    n = 1

    # --- Calculation ---
    # 1. Calculate the epsilon parameter
    epsilon = d_c - d

    # 2. Use the first-order formula for alpha
    # α = (4 - n) / (2 * (n + 8)) * ε
    numerator = 4 - n
    denominator_term_1 = 2
    denominator_term_2 = n + 8
    denominator = denominator_term_1 * denominator_term_2
    
    alpha = (numerator / denominator) * epsilon

    # --- Output ---
    print("This script calculates the specific heat scaling exponent α for d=3 using the epsilon expansion.")
    print(f"We assume the Ising model, which has n={n} component(s) for the order parameter.")
    print("-" * 30)
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The epsilon parameter is calculated as: ε = d_c - d = {d_c} - {d} = {epsilon}.")
    print("\nThe first-order expansion formula for α is: α = (4 - n) / (2 * (n + 8)) * ε")
    print("\nPlugging in the values:")
    print(f"α = ({4} - {n}) / ({denominator_term_1} * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / ({denominator_term_1} * {denominator_term_2}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    print(f"α = {numerator/denominator} * {epsilon}")
    print(f"The final result is α = {alpha}")
    print(f"(Note: This is equivalent to the fraction 1/6)")

if __name__ == '__main__':
    calculate_alpha_exponent()