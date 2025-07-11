import math
from fractions import Fraction

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau Link.
    """
    # Step 1: Define the problem data
    # Weights of the weighted projective space P(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]
    
    # The polynomial is P = z_1^8*z_3 + z_1^4*z_2^3*z_3 + z_1*z_2^7 + z_1*z_2*z_3*z_4*z_5 + z_2*z_3^4 + z_4^3*z_5 + z_5^3
    # For a Calabi-Yau hypersurface, the degree 'd' is the sum of the weights.
    # Let's confirm this. The sum of weights is:
    # 22 + 29 + 49 + 50 + 75 = 225. So, d = 225.
    # We can check the weighted degrees of the monomials:
    # deg(z_1^8*z_3) = 8*22 + 1*49 = 176 + 49 = 225
    # deg(z_1^4*z_2^3*z_3) = 4*22 + 3*29 + 1*49 = 88 + 87 + 49 = 224 (This appears to be a typo in the polynomial)
    # All other monomials have a degree of 225. We will proceed with d=225, which is consistent
    # with the Calabi-Yau condition and the majority of the polynomial's terms.
    
    d = sum(weights)

    # Step 2: Calculate the least common multiple (LCM) of the weights.
    def gcd(a, b):
        """Computes the greatest common divisor of a and b."""
        return math.gcd(a, b)

    def lcm(a, b):
        """Computes the least common multiple of a and b."""
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b) if a != 0 and b != 0 else 0

    def lcm_list(numbers):
        """Computes the least common multiple of a list of numbers."""
        result = 1
        for number in numbers:
            result = lcm(result, number)
        return result

    lcm_of_weights = lcm_list(weights)

    # Step 3: Calculate the Crawley-Nordström invariant
    # The invariant c(X) is defined as d / lcm(weights)
    invariant = Fraction(d, lcm_of_weights)
    numerator = invariant.numerator
    denominator = invariant.denominator

    # Step 4: Print the results as requested
    print("This script calculates the Crawley-Nordström invariant for the specified Calabi-Yau Link.")
    print("-" * 50)
    
    print("Step 1: Determine the degree 'd' of the polynomial.")
    print("For a Calabi-Yau hypersurface, the degree 'd' equals the sum of the weights.")
    print(f"d = {' + '.join(map(str, weights))} = {d}")
    print("-" * 50)

    print("Step 2: Calculate the least common multiple (LCM) of the weights.")
    print(f"Weights = {tuple(weights)}")
    print(f"LCM{tuple(weights)} = {lcm_of_weights}")
    print("-" * 50)

    print("Step 3: Calculate the Crawley-Nordström invariant (c).")
    print("The invariant is defined as c = d / LCM.")
    print(f"c = {d} / {lcm_of_weights}")
    print("\nSimplifying the fraction, we get the final result:")
    print(f"c = {numerator} / {denominator}")
    print("-" * 50)

if __name__ == "__main__":
    solve_crawley_nordstrom()