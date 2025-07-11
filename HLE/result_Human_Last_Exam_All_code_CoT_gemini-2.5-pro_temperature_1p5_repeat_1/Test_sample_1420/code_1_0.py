import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha (α) using the
    first-order epsilon expansion in renormalization group theory.
    """
    # 1. Define the parameters
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension for the O(N) model
    N = 1  # Number of components for the order parameter (Ising model)

    # 2. Calculate epsilon (ϵ)
    epsilon = d_c - d

    # Use sympy for symbolic representation to keep the fraction
    N_sym = sympy.Symbol('N')
    epsilon_sym = sympy.Symbol('epsilon')

    # 3. Define the formula for alpha
    alpha_formula_sym = (4 - N_sym) / (2 * (N_sym + 8)) * epsilon_sym

    # Substitute the values to get the result as a fraction
    alpha_val_sym = alpha_formula_sym.subs({N_sym: N, epsilon_sym: epsilon})
    alpha_val_float = float(alpha_val_sym)

    # 4. Print the explanation and the step-by-step calculation
    print("Calculating the specific heat scaling exponent α for d=3 using the first-order epsilon expansion.")
    print("-" * 70)
    print(f"Given parameters:")
    print(f"  - Spatial dimension (d): {d}")
    print(f"  - Upper critical dimension (d_c): {d_c}")
    print(f"  - Order parameter components (N, for Ising model): {N}\n")

    print(f"Step 1: Calculate the expansion parameter epsilon (ϵ)")
    print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}\n")

    print("Step 2: Apply the first-order expansion formula for α")
    print(f"α = (4 - N) / (2 * (N + 8)) * ϵ\n")

    print("Step 3: Substitute the values into the formula")
    print(f"α = (4 - {N}) / (2 * ({N} + 8)) * {epsilon}")
    print(f"α = ({4 - N}) / (2 * {N + 8}) * {epsilon}")
    print(f"α = {4 - N} / {2 * (N + 8)}")
    print(f"α = {alpha_val_sym}\n")

    print("Final result:")
    print(f"The scaling exponent α for d=3 is {alpha_val_sym}, which is approximately {alpha_val_float:.4f}.")

if __name__ == "__main__":
    calculate_alpha_exponent()
