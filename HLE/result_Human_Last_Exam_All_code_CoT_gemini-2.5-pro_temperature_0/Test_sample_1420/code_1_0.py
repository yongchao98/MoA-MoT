import sympy

def calculate_specific_heat_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.
    """
    # 1. Define the parameters
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension for the O(N) model
    N = 1      # Number of components for the order parameter (Ising model)

    # 2. Calculate epsilon
    epsilon = d_c - d

    # 3. Use the first-order formula for alpha
    # alpha = (4 - N) / (2 * (N + 8)) * epsilon
    # We use sympy to keep the result as a fraction for precision
    num = 4 - N
    den = 2 * (N + 8)
    alpha = sympy.Rational(num, den) * epsilon

    # 4. Print the explanation and the calculation steps
    print("Calculation of the specific heat exponent alpha using the epsilon expansion:")
    print(f"The spatial dimension is d = {d}.")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The expansion parameter is epsilon = d_c - d = {d_c} - {d} = {epsilon}.")
    print(f"Assuming the Ising model universality class, the number of components is N = {N}.")
    print("\nThe first-order formula for alpha is: (4 - N) / (2 * (N + 8)) * epsilon")
    print("\nSubstituting the values into the equation:")
    print(f"alpha = ({4} - {N}) / (2 * ({N} + {8})) * {epsilon}")
    print(f"alpha = {num} / ({2} * {N+8}) * {epsilon}")
    print(f"alpha = {num} / {den} * {epsilon}")
    print(f"alpha = {alpha}")
    print(f"As a decimal, alpha is approximately: {float(alpha):.4f}")

if __name__ == '__main__':
    calculate_specific_heat_exponent()