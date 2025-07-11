import sys

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.

    In renormalization group theory, the epsilon expansion is used to calculate
    critical exponents near the upper critical dimension (d_c). For the Ising
    universality class, d_c = 4. The expansion is done in terms of the
    small parameter epsilon = d_c - d.

    The formula for the specific heat exponent alpha to first order in epsilon
    for an n-component model is:
    alpha = (4 - n) / (2 * (n + 8)) * epsilon
    """

    # --- Parameters ---
    # Spatial dimension
    d = 3
    # Upper critical dimension for this model class
    d_c = 4
    # Number of components of the order parameter (n=1 for the Ising model)
    n = 1

    # --- Calculation ---
    # 1. Calculate the expansion parameter epsilon
    epsilon = d_c - d

    # 2. Calculate the numerator and denominator of the formula
    numerator = 4 - n
    denominator = 2 * (n + 8)

    # 3. Calculate alpha
    alpha = (numerator / denominator) * epsilon

    # --- Output ---
    print(f"Calculating the specific heat exponent α for d={d} using the ε-expansion.")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The expansion parameter is ε = d_c - d = {d_c} - {d} = {epsilon}.")
    print("\nThe first-order formula for α is: α = (4 - n) / (2 * (n + 8)) * ε")
    print(f"Assuming the Ising model (n={n}), we substitute the values:")
    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    print(f"α = {numerator/denominator} * {epsilon}")
    print(f"α = {alpha}")
    print(f"\nThe result to first order in epsilon is 1/6, which is approximately {alpha:.4f}.")
    
    # Writing the final answer to a file as per the convention
    # For this interactive session, we'll print it in the required format.
    # Note: Using sys.stdout to ensure it's captured if this is run as a script.
    sys.stdout.write(f"\n<<<{alpha:.3f}>>>")

calculate_alpha_exponent()