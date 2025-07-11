import fractions

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent α for d=3 using the
    first-order epsilon expansion for the Ising model.
    """
    # 1. Define parameters
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension
    
    # Since the model is not specified, we assume the Ising model (n=1),
    # a standard choice for such problems.
    n = 1

    # 2. Calculate the expansion parameter epsilon
    epsilon = d_c - d

    # 3. Explain the method and parameters
    print("To find the scaling exponent α for the specific heat at d=3, we use the epsilon (ε) expansion.")
    print(f"The calculation is performed near the upper critical dimension, d_c = {d_c}.")
    print(f"The expansion parameter is ε = d_c - d = {d_c} - {d} = {epsilon}.\n")

    print("The formula for α to the first order in ε for an O(n)-symmetric model is:")
    print("α = (4 - n) * ε / (2 * (n + 8))\n")

    print(f"We assume the Ising model, where the number of order parameter components is n = {n}.")
    print("Substituting n=1 and ε=1 into the formula, we get the equation:\n")

    # 4. Perform the calculation
    numerator = (4 - n) * epsilon
    denominator = 2 * (n + 8)
    
    # Represent the result as a fraction for precision
    alpha_frac = fractions.Fraction(numerator, denominator)
    alpha_float = float(alpha_frac)

    # 5. Print the step-by-step calculation as requested
    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = ({4 - n}) / (2 * ({n + 8})) * {epsilon}")
    print(f"α = {numerator} / {denominator}\n")

    print(f"The calculated value for the scaling exponent is α = {alpha_frac}.")
    print(f"As a decimal, this is approximately α ≈ {alpha_float:.4f}.")

if __name__ == '__main__':
    calculate_alpha_exponent()