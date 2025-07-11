import math

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau link.
    """
    # Step 1: Define weights and number of variables
    weights = [22, 29, 49, 50, 75]
    n = len(weights)

    # Step 2: Determine the degree 'd'
    # Based on the analysis in the text above, the degree is 225.
    d = 225

    # Step 3 & 4: Apply the formula and calculate
    # Numerator is d^(n-1)
    # Denominator is the product of weights
    
    numerator = d**(n - 1)
    
    denominator = 1
    for w in weights:
        denominator *= w

    # Simplify the fraction by dividing by the greatest common divisor (GCD)
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    # Step 5: Present the result
    print("The Crawley-Nordström invariant is calculated using the formula: d^(n-1) / (w1 * ... * wn)")
    print(f"For degree d = {d} and n = {n} variables with weights = {weights}, the equation is:")
    
    # Building the string for the denominator product
    denominator_str = " * ".join(map(str, weights))
    
    print(f"{d}^{n-1} / ({denominator_str})")
    print(f"= {numerator} / {denominator}")
    print(f"= {simplified_numerator} / {simplified_denominator}")

calculate_crawley_nordstrom_invariant()