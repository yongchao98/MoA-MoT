import sys

def calculate_alpha_exponent():
    """
    Calculates the specific heat exponent alpha using the first-order epsilon expansion.
    """
    # Parameters for the calculation
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension for the O(N) model
    N = 1      # Number of order parameter components (Ising model)

    # Calculate epsilon
    epsilon = d_c - d

    # Calculate alpha using the first-order expansion formula
    numerator = 4 - N
    denominator = 2 * (N + 8)
    alpha = float(numerator) / float(denominator) * epsilon

    # Print the explanation and the steps of the calculation
    print("Calculation of the specific heat exponent 'alpha' using the first-order epsilon expansion.")
    print("The formula is: alpha = (4 - N) / (2 * (N + 8)) * epsilon")
    print("-" * 50)
    print(f"Given values:")
    print(f"Spatial dimension, d = {d}")
    print(f"Upper critical dimension, d_c = {d_c}")
    print(f"Number of components (Ising model), N = {N}")
    print("-" * 50)
    print(f"Step 1: Calculate epsilon")
    print(f"epsilon = d_c - d = {d_c} - {d} = {epsilon}")
    print("-" * 50)
    print(f"Step 2: Substitute values into the formula for alpha")
    # Show the final equation with all numbers substituted as requested
    print(f"alpha = ({4} - {N}) / ({2} * ({N} + {8})) * {epsilon}")
    print(f"alpha = {numerator} / ({2} * {N+8}) * {epsilon}")
    print(f"alpha = {numerator} / {denominator}")
    print("-" * 50)
    print(f"The resulting value for the specific heat exponent alpha is: {alpha}")


calculate_alpha_exponent()