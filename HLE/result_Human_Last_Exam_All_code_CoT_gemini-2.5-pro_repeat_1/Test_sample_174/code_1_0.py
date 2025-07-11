import math
from fractions import Fraction

def compute_bernoulli_numbers(n):
    """
    Computes Bernoulli numbers B_0, B_1, ..., B_n using the recurrence relation.
    sum_{k=0 to m} C(m+1, k) * B_k = 0 for m >= 1.
    """
    B = [Fraction(0)] * (n + 1)
    B[0] = Fraction(1)
    for m in range(1, n + 1):
        s = Fraction(0)
        for k in range(m):
            # math.comb is available in Python 3.8+
            # For older versions, one could implement it as factorial(n)/(factorial(k)*factorial(n-k))
            s += math.comb(m + 1, k) * B[k]
        B[m] = -s / (m + 1)
    return B

def main():
    """
    Calculates the orbifold Euler characteristic of the space of smooth plane quartics.
    """
    g = 3
    print(f"The genus of a smooth plane quartic curve is g = {g}.")
    print("-" * 30)

    # 1. Compute Bernoulli number B_{2g} = B_6
    n_bernoulli = 2 * g
    bernoulli_numbers = compute_bernoulli_numbers(n_bernoulli)
    B_2g = bernoulli_numbers[n_bernoulli]
    print(f"Step 1: Compute the Bernoulli number B_{2g} = B_{n_bernoulli}.")
    print(f"B_{n_bernoulli} = {B_2g}")
    print("-" * 30)

    # 2. Compute chi_orb(M_3)
    print(f"Step 2: Compute chi_orb(M_{g}) using the Harer-Zagier formula.")
    print(f"Formula: chi_orb(M_g) = -B_{2g} / (2g * (2 - 2g))")
    
    numerator = -B_2g
    denominator = 2 * g * (2 - 2 * g)
    chi_M3 = Fraction(numerator, denominator)
    
    print(f"chi_orb(M_{g}) = -({B_2g}) / ({2*g} * ({2 - 2*g})) = -({B_2g}) / ({denominator}) = {chi_M3}")
    print("-" * 30)
    
    # 3. Use the known value for chi_orb(H_3)
    chi_H3 = Fraction(-1, 120)
    print(f"Step 3: Use the known value for the hyperelliptic locus chi_orb(H_{g}).")
    print(f"chi_orb(H_{g}) = {chi_H3}")
    print("-" * 30)
    
    # 4. Calculate the final result
    print(f"Step 4: Compute the final result.")
    print(f"chi_orb([U/G]) = chi_orb(M_{g}) - chi_orb(H_{g})")
    
    result = chi_M3 - chi_H3
    
    # Print the final equation with all the numbers
    print("\nFinal Equation:")
    print(f"{chi_M3} - ({chi_H3}) = {chi_M3} + {-chi_H3} = {result}")
    
    print("\nThe orbifold Euler characteristic of the quotient stack [U/G] is:")
    print(result)

if __name__ == "__main__":
    main()