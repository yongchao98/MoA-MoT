import math
from fractions import Fraction

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the epsilon expansion.
    """
    # --- Parameters ---
    # Spatial dimension
    d = 3
    # Upper critical dimension for the n-vector model
    d_c = 4
    # Number of components of the order parameter.
    # We assume the Ising model universality class, where n=1.
    n = 1

    # --- Calculation ---
    # 1. Calculate epsilon
    epsilon = d_c - d

    # 2. Use the first-order epsilon expansion formula for alpha
    # Formula: alpha = (4 - n) / (2 * (n + 8)) * epsilon
    numerator = 4 - n
    denominator = 2 * (n + 8)
    alpha = (numerator / denominator) * epsilon
    
    # --- Output ---
    print("This script calculates the specific heat exponent α for d=3 using the first-order epsilon expansion.")
    print("\nStep 1: Define the parameters.")
    print(f"Spatial dimension, d = {d}")
    print(f"Upper critical dimension, d_c = {d_c}")
    print(f"Epsilon, ϵ = d_c - d = {d_c} - {d} = {epsilon}")
    print(f"Assuming the Ising model, the number of components of the order parameter is n = {n}.")

    print("\nStep 2: Apply the formula for α.")
    print("The formula is: α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("\nSubstituting the values into the equation:")
    print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / (2 * {n+8}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")

    # Represent as a fraction for precision
    alpha_fraction = Fraction(int(numerator * epsilon), int(denominator))
    print(f"\nThe exact result to first order in epsilon is: {alpha_fraction.numerator}/{alpha_fraction.denominator}")
    
    print(f"\nAs a decimal, the scaling exponent α is approximately: {alpha:.3f}")

if __name__ == "__main__":
    calculate_alpha_exponent()