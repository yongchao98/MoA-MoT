import math
from fractions import Fraction

def solve_integral():
    """
    This function calculates the integral by first computing the sum S_50,
    and then using it to find the value of the integral.
    It then prints the final equation and the result.
    """
    n = 50

    # Calculate S_n = sum_{k} (-1)^k comb(n, k) comb(n, 3k-n)
    # This sum corresponds to the coefficient of z^(2n) in (z-1)^n(z^3+1)^n
    # with an extra factor of (-1)^n.
    # We want S_mine = [z^{2n}](z-1)^n(z^3+1)^n = (-1)^n * S_n.
    
    s = 0
    k_min = (n + 2) // 3
    k_max = (2 * n) // 3
    
    # Store terms of the summation to display them
    sum_terms = []
    
    for k in range(k_min, k_max + 1):
        try:
            term = math.comb(n, k) * math.comb(n, 3 * k - n)
            if k % 2 == 1:
                s -= term
                sum_terms.append(f"-C({n},{k})*C({n},{3*k-n})")
            else:
                s += term
                sum_terms.append(f"+C({n},{k})*C({n},{3*k-n})")
        except ValueError:
            # comb(n, r) raises ValueError if r < 0 or r > n
            pass
            
    # The sum we want, S_50_mine, is (-1)^n times the sum 's' calculated above.
    # As n=50 is even, (-1)^n = 1, so S_50_mine = s.
    S_50_mine = s

    # The formula for the integral I_n for n=50 is:
    # I_n = (-pi / 2^n) * S_n_mine, because i^n = i^50 = i^2 = -1.
    
    print("The integral is I = integral from 0 to pi of (sin(4x) - sin(2x))^50 dx.")
    print("Using complex analysis, this integral can be shown to be I = (-pi / 2^50) * S_50.")
    print("S_50 is a sum of binomial coefficients: S_50 = sum_{k=17 to 33} (-1)^k * C(50,k) * C(50, 3k-50).")
    print(f"Let's compute S_50:")
    print(f"S_50 = {' '.join(sum_terms)}")
    print(f"S_50 = {S_50_mine}")
    
    numerator = -S_50_mine
    denominator = 2**n
    
    final_fraction = Fraction(numerator, denominator)
    
    print("\nSo the final equation for the integral is:")
    print(f"I = (-pi / {2**n}) * ({S_50_mine})")
    print(f"I = pi * ({-S_50_mine} / {2**n})")
    print(f"I = pi * ({final_fraction.numerator} / {final_fraction.denominator})")
    
    print("\nThe value of the integral is pi times the following fraction:")
    print(f"{final_fraction.numerator}/{final_fraction.denominator}")
    
    # The final answer in the requested format
    # Note: the prompt implies the answer itself should be a fraction,
    # suggesting we should output I/pi.
    print(f"\nFinal Answer (value of I/pi):")
    print(f"<<<{final_fraction.numerator}/{final_fraction.denominator}>>>")

solve_integral()