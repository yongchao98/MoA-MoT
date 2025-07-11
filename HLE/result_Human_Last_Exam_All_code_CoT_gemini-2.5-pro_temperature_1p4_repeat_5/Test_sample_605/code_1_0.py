import math

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau hypersurface.
    """
    # The weights of the ambient projective space P^4(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]

    # The Calabi-Yau condition states that the degree d is the sum of the weights.
    d = sum(weights)

    # Numerator of the invariant formula is d^3.
    numerator = d**3

    # Denominator of the invariant formula is the product of the weights.
    denominator = 1
    for w in weights:
        denominator *= w

    # To simplify the fraction, find the greatest common divisor (GCD).
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    # Print the step-by-step calculation
    print("The formula for the Crawley-Nordström invariant is c(X) = d^3 / (w1 * w2 * w3 * w4 * w5)")
    print(f"The weights are w = {tuple(weights)}")
    print("The degree d is the sum of the weights (Calabi-Yau condition):")
    weight_sum_str = " + ".join(map(str, weights))
    print(f"d = {weight_sum_str} = {d}\n")
    
    print("Calculating the numerator and denominator:")
    print(f"Numerator = d^3 = {d}^3 = {numerator}")
    weight_prod_str = " * ".join(map(str, weights))
    print(f"Denominator = {weight_prod_str} = {denominator}\n")
    
    print("The invariant as a fraction is:")
    print(f"c(X) = {numerator} / {denominator}\n")
    
    print("Simplifying the fraction gives the final result:")
    print(f"c(X) = {simplified_numerator} / {simplified_denominator}")

calculate_crawley_nordstrom_invariant()
<<<6075/62524>>>