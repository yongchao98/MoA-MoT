import math

def solve_pm(m):
    """
    Calculates the probability Pm for a given positive integer m.

    Args:
      m: A positive integer.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.")
        return

    # Calculate the components of the numerator of Pm
    num_factor1 = m + 1
    num_factor2 = m + 2
    numerator = num_factor1 * num_factor2

    # Calculate the components of the denominator of Pm
    den_factor1 = 2
    den_factor2 = 2 * m + 1
    den_factor3 = 4 * m + 1
    denominator = den_factor1 * den_factor2 * den_factor3

    # To present the result as a simplified fraction
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    # Output the results
    print(f"For m = {m}:")
    print("The formula for the probability Pm is (m+1)(m+2) / [2(2m+1)(4m+1)]")
    print("\nCalculation:")
    print(f"Numerator = ({m}+1) * ({m}+2) = {num_factor1} * {num_factor2} = {numerator}")
    print(f"Denominator = 2 * (2*{m}+1) * (4*{m}+1) = {den_factor1} * {den_factor2} * {den_factor3} = {denominator}")
    print(f"\nPm = {numerator}/{denominator}")
    print(f"The simplified probability is: {simplified_num}/{simplified_den}")


# Example usage with m=1
# You can change this value to test other positive integers for m.
m_value = 1
solve_pm(m_value)

# Example usage with m=2
m_value = 2
solve_pm(m_value)
