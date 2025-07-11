import math
from itertools import combinations
from fractions import Fraction

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau hypersurface
    in a weighted projective space.
    """
    # 1. Define the weights and verify the Calabi-Yau condition.
    weights = [22, 29, 49, 50, 75]
    
    # For a Calabi-Yau hypersurface, the degree 'd' equals the sum of the weights.
    d = sum(weights)

    print("Step 1: Define the parameters from the problem.")
    print(f"The weights of the ambient space are w = {weights}")
    print(f"The degree of the Calabi-Yau hypersurface is d = sum(w) = {d}\n")
    
    # 2. Calculate the necessary components for the formula.
    # Sigma_2 is the sum of pairwise products of weights.
    sigma_2 = sum(p[0] * p[1] for p in combinations(weights, 2))
    
    # The product of all weights.
    prod_weights = math.prod(weights)
    
    # 3. Apply the formula for the Crawley-Nordstrom invariant.
    # CN = (d * sigma_2) / (product of weights)
    numerator = d * sigma_2
    denominator = prod_weights

    # Simplify the fraction for the final answer.
    final_fraction = Fraction(numerator, denominator)

    print("Step 2: Use the formula for the Crawley-Nordstrom invariant:")
    print("CN = (d * Σ₂) / (w1*w2*w3*w4*w5)")
    print("Where Σ₂ is the sum of pairwise products of weights.\n")

    print("Step 3: Calculate each part of the equation.")
    print(f"  d = {d}")
    print(f"  Σ₂ = {sigma_2}")
    print(f"  w1*w2*w3*w4*w5 = {denominator}\n")

    print("Step 4: Substitute the values into the formula.")
    print(f"  CN = ({d} * {sigma_2}) / {denominator}")
    print(f"     = {numerator} / {denominator}\n")
    
    print("Step 5: Simplify the fraction to get the final answer.")
    print(f"  Crawley-Nordstrom Invariant = {final_fraction.numerator} / {final_fraction.denominator}")

solve_crawley_nordstrom()