import sympy

def calculate_specific_heat_exponent():
    """
    Calculates the specific heat scaling exponent α using the first-order
    epsilon expansion for the 3D Ising universality class.
    """
    # Parameters for the epsilon expansion
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension
    n = 1  # Number of order parameter components (n=1 for the Ising model)

    # Calculate epsilon
    epsilon = d_c - d

    # The formula for alpha to the first order in epsilon is:
    # α = (4 - n) / (2 * (n + 8)) * ϵ
    numerator_val = 4 - n
    denominator_val = 2 * (n + 8)

    # Calculate alpha as a fraction for precision
    alpha_frac = sympy.Rational(numerator_val, denominator_val) * epsilon
    alpha_float = float(alpha_frac)

    # Print the explanation and the step-by-step calculation
    print("Calculating the specific heat exponent α for d=3 using the ϵ-expansion.")
    print("We assume the system is in the Ising universality class (order parameter components n=1).")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The spatial dimension is d = {d}.")
    print(f"Thus, ϵ = d_c - d = {d_c} - {d} = {epsilon}.")
    print("\nThe first-order expansion formula for α is: α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("\nSubstituting the values n=1 and ϵ=1:")
    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator_val} / (2 * {n+8}) * {epsilon}")
    print(f"α = {numerator_val} / {denominator_val} * {epsilon}")
    print(f"α = {alpha_frac}")
    print(f"\nAs a decimal, the value is approximately: {alpha_float:.4f}")

if __name__ == '__main__':
    calculate_specific_heat_exponent()
    # The final answer in fraction form is 1/6
    final_answer = 1/6
    # To conform to the output format, we can output the float value
    # print(f'<<<{final_answer:.4f}>>>')