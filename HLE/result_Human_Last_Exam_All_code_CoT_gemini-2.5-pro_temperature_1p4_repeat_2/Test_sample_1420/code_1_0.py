import fractions

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha (α) using the epsilon expansion.

    In the context of renormalization group theory, the epsilon expansion is a technique
    used to calculate critical exponents near the upper critical dimension (d_c).
    For the n-vector model, the upper critical dimension is d_c = 4.

    The scaling exponent α for the specific heat is given to the first order in epsilon (ϵ) by the formula:
    α = (4 - n) / (2 * (n + 8)) * ϵ

    We will solve this for a system with spatial dimension d=3. The canonical example
    for this is the Ising model, where the order parameter has n=1 component.
    """
    # Parameters
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension
    n = 1  # Number of components of the order parameter (for the Ising model)

    # Calculate epsilon
    epsilon = d_c - d

    # Numerator and denominator of the main fraction
    numerator = 4 - n
    denominator = 2 * (n + 8)

    # Calculate alpha as a fraction for precision
    alpha_fraction = fractions.Fraction(numerator * epsilon, denominator)

    # Calculate alpha as a float for a numerical approximation
    alpha_float = float(alpha_fraction)

    print("Calculation of the Specific Heat Scaling Exponent α")
    print("--------------------------------------------------")
    print(f"Spatial dimension, d = {d}")
    print(f"Upper critical dimension, d_c = {d_c}")
    print(f"Order parameter components (Ising model), n = {n}")
    print(f"Epsilon, ϵ = d_c - d = {d_c} - {d} = {epsilon}\n")

    print("The first-order formula for α is: α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("Plugging in the values, the equation becomes:")
    # Printing the equation with all the numbers
    print(f"α = ({4} - {n}) / ({2} * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / {denominator}")
    print(f"α = {alpha_fraction}\n")

    print("Final result:")
    print(f"The scaling exponent α is approximately: {alpha_float:.5f}")

calculate_alpha_exponent()