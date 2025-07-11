import math

def find_Pm(m):
    """
    Calculates the probability Pm for a given integer m.

    The formula is derived from combinatorial analysis:
    Number of successful pairs (i,j): N = C(m+2, 2)
    Total number of pairs (i,j): T = C(4m+2, 2)
    Pm = N / T
    """
    if not isinstance(m, int) or m <= 0:
        print("m must be a positive integer.")
        return

    # Numerator N = (m+1)(m+2)/2
    num_val_m_plus_1 = m + 1
    num_val_m_plus_2 = m + 2
    numerator = (num_val_m_plus_1 * num_val_m_plus_2) // 2
    
    # Denominator T = (2m+1)(4m+1)
    den_val_2m_plus_1 = 2 * m + 1
    den_val_4m_plus_1 = 4 * m + 1
    denominator = den_val_2m_plus_1 * den_val_4m_plus_1

    # Print the equation with substituted values
    print(f"For m = {m}:")
    print(f"Number of valid pairs (N) = (({m}+1) * ({m}+2)) / 2 = ({num_val_m_plus_1} * {num_val_m_plus_2}) / 2 = {numerator}")
    print(f"Total number of pairs (T) = (2*{m}+1) * (4*{m}+1) = {den_val_2m_plus_1} * {den_val_4m_plus_1} = {denominator}")
    
    # Simplify the fraction
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print(f"The probability P_{m} = N / T = {numerator} / {denominator}")
    print(f"Simplified P_{m} = {simplified_num} / {simplified_den}")
    print(f"As a decimal, P_{m} approx {numerator/denominator:.6f}")


# Example usage, you can change the value of m
m_value = 2
find_Pm(m_value)

m_value = 1
find_Pm(m_value)