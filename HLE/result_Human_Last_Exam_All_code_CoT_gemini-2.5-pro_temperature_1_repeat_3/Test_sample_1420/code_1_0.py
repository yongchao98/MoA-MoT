import sys

def calculate_alpha_exponent():
    """
    Calculates the specific heat exponent alpha (α) using the first-order epsilon expansion.
    """
    # Parameters for the epsilon expansion
    d_c = 4  # Upper critical dimension
    d = 3    # Spatial dimension
    n = 1    # Number of order parameter components (for the Ising model)

    # Step 1: Calculate epsilon
    epsilon = d_c - d

    # Step 2: Calculate the components of the formula for alpha
    numerator = 4 - n
    inner_denominator = n + 8
    denominator = 2 * inner_denominator

    # Step 3: Calculate alpha to first order in epsilon
    alpha = (numerator / denominator) * epsilon

    # Print the explanation and the step-by-step calculation
    print("Calculating the specific heat exponent alpha (α) using the first-order epsilon expansion.")
    print("-" * 70)
    print(f"This calculation assumes the Ising model universality class (n=1).")
    print(f"Upper critical dimension (d_c): {d_c}")
    print(f"Spatial dimension (d): {d}")
    print(f"Number of components (n): {n}")
    print(f"Epsilon (ε = d_c - d): {d_c} - {d} = {epsilon}")
    print("\nFormula for alpha: α = (4 - n) / (2 * (n + 8)) * ε")
    
    # Print the equation with all numbers substituted
    print("\nSubstituting the values into the formula:")
    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / (2 * {inner_denominator}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    
    # Print the final result
    print(f"\nThe calculated scaling exponent alpha is: {alpha}")
    
    # Required for the final answer format
    # Redirect final answer to a specific stream if needed, or just print
    # For this problem, we will print it at the end to be captured.
    print(f"<<<{alpha}>>>", file=sys.stdout)

if __name__ == '__main__':
    calculate_alpha_exponent()