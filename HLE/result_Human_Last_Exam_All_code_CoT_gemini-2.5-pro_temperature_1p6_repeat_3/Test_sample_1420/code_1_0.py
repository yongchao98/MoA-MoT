import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the first-order
    epsilon expansion for the 3D Ising model.
    """
    # Define parameters
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension for the Ising model
    n = 1  # Number of components of the order parameter (n=1 for Ising)

    # Calculate epsilon
    epsilon = d_c - d

    # Numerator and denominator components of the formula for alpha
    num_n_term = 4 - n
    den_n_term = n + 8
    den_factor = 2

    # Perform the calculation symbolically to get a fraction
    alpha_frac = sympy.Rational(num_n_term, den_factor * den_n_term) * epsilon

    # Perform the calculation numerically
    alpha_float = float(alpha_frac)

    # Print the explanation and step-by-step calculation
    print("To find the specific heat exponent α in d=3 using the epsilon expansion, we use the first-order formula for an O(n) model:")
    print("α = (4 - n) / (2 * (n + 8)) * ε")
    print("\nFor the 3D Ising model, the parameters are:")
    print(f"  - Spatial dimension d = {d}")
    print(f"  - Upper critical dimension d_c = {d_c}")
    print(f"  - Order parameter components n = {n}")
    print("\nFirst, we calculate ε:")
    print(f"  ε = d_c - d = {d_c} - {d} = {epsilon}")
    print("\nNow, we substitute n=1 and ε=1 into the formula:")
    # Using format to show each number in the equation
    print(f"α = ({4} - {n}) / ({den_factor} * ({n} + {den_n_term})) * {epsilon}")
    print(f"α = {num_n_term} / ({den_factor} * {den_n_term}) * {epsilon}")
    print(f"α = {num_n_term} / {den_factor * den_n_term}")
    print(f"\nThe result is:")
    print(f"α = {alpha_frac}")
    print(f"As a decimal, α ≈ {alpha_float}")

if __name__ == "__main__":
    calculate_alpha_exponent()