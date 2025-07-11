import math

def solve_pm():
    """
    Calculates the probability P_m for a given positive integer m.
    """
    try:
        m_str = input("Enter a positive integer m: ")
        m = int(m_str)
        if m <= 0:
            print("Error: m must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a positive integer.")
        return

    # Calculate the numerator of the probability formula P_m
    # Numerator = (m+1) * (m+2)
    num_factor1 = m + 1
    num_factor2 = m + 2
    numerator = num_factor1 * num_factor2

    # Calculate the denominator of the probability formula P_m
    # Denominator = 2 * (2m+1) * (4m+1)
    den_factor1 = 2 * m + 1
    den_factor2 = 4 * m + 1
    denominator = 2 * den_factor1 * den_factor2

    # To present the fraction in its simplest form, we find the greatest common divisor (GCD).
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor
    
    # Output the final equation and the result
    print(f"\nThe formula for P_m is (m+1)(m+2) / [2(2m+1)(4m+1)]")
    print(f"\nFor m = {m}:")
    print(f"P_{m} = ({m}+1) * ({m}+2) / (2 * (2*{m}+1) * (4*{m}+1))")
    print(f"P_{m} = {num_factor1} * {num_factor2} / (2 * {den_factor1} * {den_factor2})")
    print(f"P_{m} = {numerator} / {denominator}")
    if common_divisor > 1:
        print(f"P_{m} = {simplified_num} / {simplified_den}")
    
    probability_value = numerator / denominator
    print(f"P_{m} ≈ {probability_value}")


if __name__ == "__main__":
    solve_pm()
