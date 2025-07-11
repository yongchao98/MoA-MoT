import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.
    """
    # Define parameters as rational numbers for precise symbolic calculation
    d = sympy.Rational(3)   # Spatial dimension
    d_c = sympy.Rational(4) # Upper critical dimension for the Ising model class
    n = sympy.Rational(1)   # Number of components of the order parameter (Ising model)

    # Calculate epsilon
    epsilon = d_c - d

    # Define the formula for alpha to first order in epsilon
    # alpha = (4 - n) / (2 * (n + 8)) * epsilon
    numerator = 4 - n
    denominator = 2 * (n + 8)
    alpha = (numerator / denominator) * epsilon

    # Print the explanation and the step-by-step calculation
    print("This script calculates the specific heat scaling exponent α for the Ising model")
    print("in d=3 dimensions using the first-order epsilon expansion.")
    print("-" * 60)
    print(f"The spatial dimension is d = {int(d)}.")
    print(f"The upper critical dimension is d_c = {int(d_c)}.")
    print(f"The number of order parameter components (Ising model) is n = {int(n)}.")
    print(f"\nThe expansion parameter is calculated as ε = d_c - d = {int(d_c)} - {int(d)} = {int(epsilon)}.")
    print("\nThe first-order formula for the exponent α is: (4 - n) / (2 * (n + 8)) * ε")
    print("\nSubstituting the values into the formula:")
    # Using sympy objects to show the fractional form
    final_equation = f"α = ({int(4)} - {int(n)}) / (2 * ({int(n)} + {int(8)})) * {int(epsilon)}"
    print(final_equation)
    
    step2_equation = f"α = {int(numerator)} / (2 * {int(n+8)}) * {int(epsilon)}"
    print(step2_equation)

    step3_equation = f"α = {int(numerator)} / {int(denominator)} * {int(epsilon)}"
    print(step3_equation)
    
    final_result_fraction = f"α = {alpha}"
    print(final_result_fraction)

    final_result_decimal = f"α ≈ {float(alpha):.4f}"
    print(final_result_decimal)


if __name__ == "__main__":
    calculate_alpha_exponent()