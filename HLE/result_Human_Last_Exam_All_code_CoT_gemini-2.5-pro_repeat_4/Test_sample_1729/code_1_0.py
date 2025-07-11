import math

def find_probability_Pm(m):
    """
    Calculates and prints the probability P_m for a given positive integer m.

    The probability P_m is that a sequence a_1, ..., a_{4m+2} is an
    (i,j)-divisible sequence. This has been derived to be:
    P_m = (m+1)(m+2) / (2 * (2m+1) * (4m+1))
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.")
        return

    # Calculate the number of favorable pairs
    favorable_pairs_num = (m + 1) * (m + 2)
    favorable_pairs_den = 2
    favorable_pairs = favorable_pairs_num / favorable_pairs_den

    # Calculate the total number of pairs
    total_pairs = (2 * m + 1) * (4 * m + 1)

    print(f"For m = {m}:")
    print(f"1. The number of favorable pairs (i, j) is calculated as ((m+1)*(m+2))/2:")
    print(f"   (({m}+1)*({m}+2))/2 = ({m+1}*{m+2})/2 = {favorable_pairs_num}/2 = {int(favorable_pairs)}")
    
    print(f"2. The total number of pairs (i, j) is calculated as (2m+1)*(4m+1):")
    print(f"   ({2*m}+1)*({4*m}+1) = {2*m+1}*{4*m+1} = {total_pairs}")
    
    # Final probability calculation
    numerator = (m + 1) * (m + 2)
    denominator = 2 * (2 * m + 1) * (4 * m + 1)
    
    print("\n3. The probability P_m is the ratio of favorable pairs to total pairs.")
    print("   The final equation is P_m = (m+1)(m+2) / (2 * (2m+1) * (4m+1))")
    print(f"   Substituting m = {m}:")
    print(f"   P_{m} = (({m} + 1) * ({m} + 2)) / (2 * ({2*m} + 1) * ({4*m} + 1))")
    print(f"   P_{m} = ({m+1} * {m+2}) / (2 * {2*m+1} * {4*m+1})")
    print(f"   P_{m} = {numerator} / {denominator}")

    # Simplify the fraction
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print(f"\n   Simplified fraction: P_{m} = {simplified_num} / {simplified_den}")


# --- Execution ---
# You can change the value of m here to test with other positive integers.
m_value = 2
find_probability_Pm(m_value)