import sympy

def calculate_specific_heat_exponent():
    """
    Calculates the specific heat exponent alpha using the first-order epsilon expansion
    for the Ising model in d=3 dimensions.
    """
    # Define the parameters for the calculation
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension for the O(n) model
    n = 1  # Number of order parameter components (Ising model)

    print(f"Calculating the specific heat exponent (alpha) for a d={d} dimensional system.")
    print("We use the epsilon expansion near the upper critical dimension, d_c.")
    print(f"For the O(n) class of models, the upper critical dimension is d_c = {d_c}.")

    # Calculate epsilon
    epsilon = d_c - d
    print("\nThe expansion parameter epsilon is calculated as: epsilon = d_c - d")
    print(f"epsilon = {d_c} - {d} = {epsilon}")

    print("\nAssuming the Ising model, the number of components of the order parameter is n = 1.")
    print("The first-order epsilon expansion for alpha is given by the formula:")
    print("alpha = (4 - n) * epsilon / (2 * (n + 8))")

    # Substitute the values into the formula
    numerator = 4 - n
    denominator_part1 = n + 8
    denominator = 2 * denominator_part1
    
    # Calculate the symbolic fraction for a clear representation
    alpha_fraction = sympy.Rational(numerator * epsilon, denominator)
    alpha_float = float(alpha_fraction)

    print("\nSubstituting the values n=1 and epsilon=1:")
    # Print the equation with numbers plugged in
    print(f"alpha = ({4} - {n}) * {epsilon} / ({2} * ({n} + {8}))")
    print(f"alpha = {numerator} * {epsilon} / ({2} * {denominator_part1})")
    print(f"alpha = {numerator * epsilon} / {denominator}")
    print(f"alpha = {alpha_fraction} ~= {alpha_float:.4f}")

    print("\nThus, the scaling exponent alpha to first order in epsilon is:")
    print(f"{alpha_fraction}")

calculate_specific_heat_exponent()

alpha_val = sympy.Rational(1, 6)
# Final answer block as requested
# We need to output the numerical value. A fraction is the exact value from the first-order expansion.
print(f'<<<{alpha_val}>>>')
