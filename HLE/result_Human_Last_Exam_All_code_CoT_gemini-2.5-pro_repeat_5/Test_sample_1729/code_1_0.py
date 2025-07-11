import math

def calculate_Pm(m):
    """
    Calculates the probability P_m for a given positive integer m.

    The formula for P_m is derived as (m+1)(m+2) / (2 * (2m+1) * (4m+1)).
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.")
        return

    # Components of the numerator
    m_plus_1 = m + 1
    m_plus_2 = m + 2
    
    # Numerator of the final probability fraction
    numerator = m_plus_1 * m_plus_2

    # Components of the denominator
    two_m_plus_1 = 2 * m + 1
    four_m_plus_1 = 4 * m + 1

    # Denominator of the final probability fraction
    denominator = 2 * two_m_plus_1 * four_m_plus_1

    # To simplify the fraction, we find the greatest common divisor
    common_divisor = math.gcd(numerator, denominator)
    
    # Printing the derivation with the numbers for the given m
    print(f"For m = {m}:")
    print(f"The number of valid (i,j) pairs is C(m+2, 2) = C({m}+2, 2) = C({m_plus_2}, 2) = {int(numerator/2 * 2)}/2 = {int(numerator/2)}.")
    print(f"The total number of (i,j) pairs is C(4m+2, 2) = C({4*m}+2, 2) = C({4*m+2}, 2) = {int(denominator/2)}.")
    print("\nThe probability P_m is given by the formula:")
    print("      (m+1) * (m+2)")
    print("P_m = -------------------")
    print("    2 * (2m+1) * (4m+1)")
    print("\nPlugging in the numbers for the final equation:")
    print(f"      ({m}+1) * ({m}+2)      {m_plus_1} * {m_plus_2}      {numerator}")
    print(f"P_m = ----------------------- = -------------------- = -----------")
    print(f"    2 * (2*{m}+1) * (4*{m}+1)   2 * {two_m_plus_1} * {four_m_plus_1}     {denominator}")

    print(f"\nThe simplified probability is {numerator // common_divisor} / {denominator // common_divisor}.")

# Example usage for m=1
calculate_Pm(1)
