import math

def solve_series():
    """
    This function calculates the value of the series sum_{n=2 to infinity} (-1)^n * S_n / n
    where S_n is the n-th harmonic number.

    The analytical solution for the sum is: 1 - (pi^2 / 12) + ((ln(2))^2 / 2).
    This code calculates the value of this expression and prints the breakdown
    of the calculation as per the derivation.
    """

    # --- Step 1: Define constants and calculate components ---
    pi = math.pi
    ln2 = math.log(2)

    pi_squared = pi**2
    ln2_squared = ln2**2

    term1 = 1.0
    term2 = pi_squared / 12
    term3 = ln2_squared / 2

    # --- Step 2: Calculate the final result ---
    final_value = term1 - term2 + term3

    # --- Step 3: Print the explanation and results ---
    print("The exact value of the sum is given by the expression:")
    print("S = 1 - (pi^2 / 12) + (ln(2)^2 / 2)")
    print("\nLet's break down the calculation:")
    
    print(f"\n1. The first number in the equation is: {term1}")

    print(f"\n2. For the second term, - (pi^2 / 12):")
    print(f"   - pi is approximately: {pi}")
    print(f"   - pi^2 (pi squared) is: {pi_squared}")
    print(f"   - The denominator is: 12")
    print(f"   - So, (pi^2 / 12) is: {term2}")
    
    print(f"\n3. For the third term, + (ln(2)^2 / 2):")
    print(f"   - ln(2) (natural log of 2) is: {ln2}")
    print(f"   - ln(2)^2 (ln(2) squared) is: {ln2_squared}")
    print(f"   - The denominator is: 2")
    print(f"   - So, (ln(2)^2 / 2) is: {term3}")

    print("\nPutting it all together:")
    print(f"S = {term1} - {term2} + {term3}")
    print(f"The final value is: {final_value}")

if __name__ == '__main__':
    solve_series()
