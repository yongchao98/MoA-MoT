from fractions import Fraction

def calculate_alpha_exponent():
    """
    Calculates the specific heat exponent alpha using the first-order epsilon expansion.
    """
    # 1. Define the parameters
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension for the O(n) model
    n = 1      # Number of spin components (n=1 for the Ising universality class)

    print("Calculating the specific heat exponent α for a d=3 system.")
    print("This will be done using the epsilon (ε) expansion near the upper critical dimension.")
    print("-" * 50)
    print(f"The spatial dimension is d = {d}.")
    print(f"The upper critical dimension for the O(n) model is d_c = {d_c}.")
    print(f"We will assume the system is in the Ising universality class, so n = {n}.")
    print("-" * 50)

    # 2. Calculate epsilon
    epsilon = d_c - d
    print("The expansion parameter ε is defined as d_c - d.")
    print(f"ε = {d_c} - {d} = {epsilon}")
    print("-" * 50)

    # 3. State the formula and perform the calculation
    print("The first-order ε-expansion for the exponent α is:")
    print("α = (4 - n) / (2 * (n + 8)) * ε")
    print("\nSubstituting the values n=1 and ε=1:")

    # 4. Compute the result using the fractions module for an exact result
    numerator = 4 - n
    denominator = 2 * (n + 8)
    alpha_frac = Fraction(numerator, denominator) * epsilon

    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = {numerator} / (2 * {n + 8}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    
    # 5. Present the final output
    print("\nThe result is:")
    print(f"α = {alpha_frac.numerator}/{alpha_frac.denominator}")

    alpha_float = float(alpha_frac)
    print("\nAs a decimal, this is approximately:")
    print(f"α ≈ {alpha_float}")

if __name__ == '__main__':
    calculate_alpha_exponent()
    # The final numerical answer for the tag below
    final_answer = float(Fraction((4 - 1), (2 * (1 + 8))) * 1)