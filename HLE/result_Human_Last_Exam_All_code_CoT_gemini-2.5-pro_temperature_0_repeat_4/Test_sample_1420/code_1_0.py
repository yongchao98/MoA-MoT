import fractions

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.
    """
    # Parameters for the calculation
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension for the Ising universality class
    n = 1      # Number of components of the order parameter for the Ising model

    # 1. Calculate the expansion parameter epsilon
    epsilon = d_c - d

    # 2. Define the components of the formula for alpha
    # Formula: alpha = (4 - n) / (2 * (n + 8)) * epsilon
    numerator = 4 - n
    denominator_part1 = 2
    denominator_part2 = n + 8
    denominator = denominator_part1 * denominator_part2

    # 3. Calculate alpha
    alpha_float = (numerator / denominator) * epsilon
    alpha_fraction = fractions.Fraction(numerator, denominator) * epsilon

    # 4. Print the step-by-step calculation
    print("Calculating the specific heat exponent alpha using the first-order epsilon expansion.")
    print("-" * 70)
    print(f"The system's spatial dimension is d = {d}.")
    print(f"The upper critical dimension for the Ising universality class is d_c = {d_c}.")
    print(f"The number of order parameter components for the Ising model is n = {n}.")
    print(f"\nThe expansion parameter is calculated as: epsilon = d_c - d = {d_c} - {d} = {epsilon}")
    print("\nThe first-order expansion formula for alpha is: (4 - n) / (2 * (n + 8)) * epsilon")
    print("\nSubstituting the values into the formula:")
    # Showing the equation with all numbers substituted
    print(f"alpha = ({4} - {n}) / ({denominator_part1} * ({n} + {8})) * {epsilon}")
    print(f"alpha = {numerator} / ({denominator_part1} * {denominator_part2}) * {epsilon}")
    print(f"alpha = {numerator} / {denominator} * {epsilon}")
    print(f"\nThe calculated value for alpha is {alpha_float:.4f}, or exactly {alpha_fraction}.")

if __name__ == "__main__":
    calculate_alpha_exponent()
