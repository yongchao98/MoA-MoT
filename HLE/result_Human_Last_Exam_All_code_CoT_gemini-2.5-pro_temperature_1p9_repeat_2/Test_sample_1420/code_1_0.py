import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat exponent alpha using the epsilon expansion.

    In renormalization group theory, the epsilon expansion is a technique to
    calculate critical exponents as a power series in the parameter epsilon.
    """
    
    # Step 1: Define the parameters
    d_c = 4  # Upper critical dimension for the n-vector model
    d = 3    # Spatial dimension as given in the problem
    n = 1    # Number of components of the order parameter for the Ising model

    # Step 2: Calculate epsilon
    epsilon = d_c - d
    
    # Step 3: Calculate the alpha exponent using the first-order formula
    # alpha = (4 - n) * epsilon / (2 * (n + 8))
    # We will use SymPy to represent the result as an exact fraction.
    numerator = (4 - n) * epsilon
    denominator = 2 * (n + 8)
    
    alpha = sympy.Rational(numerator, denominator)

    # Step 4: Print the results and the calculation steps
    print("Calculation of the Specific Heat Exponent (α) via Epsilon Expansion")
    print("-" * 65)
    print(f"Upper critical dimension, d_c = {d_c}")
    print(f"Spatial dimension, d = {d}")
    print(f"Order parameter components (Ising model), n = {n}\n")
    
    print("The expansion parameter, epsilon (ε), is calculated as:")
    print(f"ε = d_c - d = {d_c} - {d} = {epsilon}\n")
    
    print("The first-order epsilon expansion formula for α is:")
    print("α = (4 - n) * ε / (2 * (n + 8))\n")

    print("Substituting the values into the formula:")
    final_equation = f"α = (4 - {n}) * {epsilon} / (2 * ({n} + {8}))"
    calculation_steps = f" = {4-n} / (2 * {n+8}) = {numerator} / {denominator}"
    final_result = f" = {alpha}"
    
    print(final_equation + calculation_steps + final_result)

calculate_alpha_exponent()