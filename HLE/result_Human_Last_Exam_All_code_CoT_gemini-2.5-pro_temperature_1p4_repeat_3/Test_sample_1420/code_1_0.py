import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.

    This calculation is for a system with spatial dimension d=3, near the
    upper critical dimension d_c=4. It assumes the Ising universality class,
    where the order parameter has n=1 component. The result is the first-order
    term in the epsilon expansion.
    """
    # Spatial dimension
    d = 3
    # Upper critical dimension
    d_c = 4
    # Number of components of the order parameter (n=1 for Ising model)
    n = 1

    # Calculate epsilon
    epsilon = d_c - d

    # Numerator and denominator of the formula for alpha
    numerator = 4 - n
    denominator = 2 * (n + 8)

    # Calculate alpha as a fraction and as a float
    alpha_fraction = sympy.Rational(numerator * epsilon, denominator)
    alpha_float = float(alpha_fraction)

    # Print the explanation and step-by-step calculation
    print("Calculating the specific heat scaling exponent α using the ϵ-expansion.")
    print(f"The system has spatial dimension d = {d}.")
    print(f"The upper critical dimension for this universality class is d_c = {d_c}.")
    print(f"Assuming the Ising model, the number of order parameter components is n = {n}.\n")

    print(f"First, we calculate ϵ = d_c - d:")
    print(f"ϵ = {d_c} - {d} = {epsilon}\n")

    print("Next, we use the first-order formula for α:")
    print("α = (4 - n) / (2 * (n + 8)) * ϵ\n")

    print("Substituting the values n=1 and ϵ=1:")
    print(f"α = ({4} - {n}) / ({2} * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / ({2} * {n + 8}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}\n")
    
    print("The final result is:")
    print(f"α = {alpha_fraction} ≈ {alpha_float}")

calculate_alpha_exponent()
<<<0.16666666666666666>>>