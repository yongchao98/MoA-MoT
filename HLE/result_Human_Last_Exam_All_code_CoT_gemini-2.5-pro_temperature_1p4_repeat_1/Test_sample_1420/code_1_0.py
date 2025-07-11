import sys
from fractions import Fraction

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.
    """
    # Step 1: Define the parameters for the calculation.
    # d_c is the upper critical dimension for O(n) models.
    d_c = 4
    # d is the spatial dimension given in the problem.
    d = 3
    # n is the number of components of the order parameter.
    # We assume the Ising model, where n=1.
    n = 1

    # Step 2: Calculate the expansion parameter, epsilon.
    epsilon = d_c - d

    # Step 3: Apply the first-order epsilon expansion formula for alpha.
    # The formula is: alpha = (4 - n) / (2 * (n + 8)) * epsilon
    numerator_term = 4 - n
    denominator_term = 2 * (n + 8)
    
    alpha_value = (numerator_term / denominator_term) * epsilon
    
    # For a more precise representation, calculate using fractions
    alpha_fraction = Fraction(numerator_term, denominator_term) * epsilon

    # Step 4: Print the process and the result.
    print("Calculation of the Specific Heat Scaling Exponent α")
    print("-" * 50)
    print(f"This calculation uses the first-order epsilon (ϵ) expansion for the Ising model (n=1).")
    print(f"Upper Critical Dimension (d_c): {d_c}")
    print(f"System Spatial Dimension (d): {d}")
    print(f"Order Parameter Components (n): {n}")
    print("\nFirst, we calculate ϵ:")
    print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}")
    
    print("\nNext, we use the formula for α:")
    print("α = (4 - n) / (2 * (n + 8)) * ϵ")
    
    print("\nSubstituting the values into the equation:")
    # The final equation with each number printed out
    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator_term} / {denominator_term} * {epsilon}")
    
    print(f"\nThe result as a fraction is: {alpha_fraction.numerator}/{alpha_fraction.denominator}")
    print(f"The numerical value for α is: {alpha_value}")
    
    # Return the final numerical answer in the required format
    # Redirecting to stdout to be captured
    sys.stdout.write(f"\n<<<{alpha_value}>>>")

calculate_alpha_exponent()