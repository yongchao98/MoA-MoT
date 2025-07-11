import sympy

def calculate_specific_heat_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the first-order
    epsilon expansion for the Ising universality class in d=3.
    """
    # Define the parameters for the calculation
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension for the phi^4 model
    n = 1      # Number of components for the order parameter (n=1 for the Ising model)

    # 1. Calculate the expansion parameter, epsilon
    epsilon = d_c - d

    # 2. Use the first-order epsilon expansion formula for alpha
    # Formula: alpha = ( (4 - n) / (2 * (n + 8)) ) * epsilon
    numerator = 4 - n
    denominator_val = n + 8
    denominator = 2 * denominator_val
    
    # Calculate alpha using floating point numbers for the final result
    alpha_val = (numerator / denominator) * epsilon

    # Use sympy for a clean fractional representation
    alpha_frac = sympy.Rational(numerator, denominator) * epsilon

    # 3. Print the step-by-step calculation
    print("Calculation of the Specific Heat Exponent (α) via Epsilon Expansion")
    print("=" * 65)
    print(f"This calculation uses the first-order (one-loop) epsilon expansion.")
    print(f"System Parameters:")
    print(f"  - Spatial Dimension (d): {d}")
    print(f"  - Upper Critical Dimension (d_c): {d_c}")
    print(f"  - Order Parameter Components (n): {n} (for the Ising universality class)\n")

    print("Step 1: Calculate the expansion parameter ϵ")
    print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}\n")

    print("Step 2: Apply the formula for α")
    print("The general first-order formula is: α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("Substituting the values into the equation:")
    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = {numerator} / (2 * {denominator_val}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}\n")

    print("Step 3: Final Result")
    print(f"The result as a fraction is: α = {alpha_frac}")
    print(f"The result as a decimal is: α ≈ {alpha_val:.4f}")
    print("=" * 65)

calculate_specific_heat_exponent()