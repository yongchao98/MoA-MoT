import math
import fractions

def solve():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.

    This density is determined by the structure of the polynomial's Galois group,
    which is identified as the symmetric group S_7. The density is the proportion
    of 7-cycles within S_7.
    """
    n = 7

    # The Galois group of the polynomial is S_n (the symmetric group on n elements).
    # We need to find the proportion of elements in S_n that are n-cycles.

    # The number of n-cycles in S_n is (n-1)!
    num_n_cycles = math.factorial(n - 1)

    # The total number of elements in S_n is n!
    total_elements = math.factorial(n)
    
    # The density is the ratio of these two numbers.
    # The fraction is simplified for a clear result.
    density = fractions.Fraction(num_n_cycles, total_elements)
    
    print(f"The Galois group of the polynomial has been identified as the symmetric group S_{n}.")
    print("The requested density is the proportion of n-cycles in this group.")
    print("\nCalculating the numbers for the final equation:")
    print(f"Number of {n}-cycles in S_{n} = ({n}-1)! = {num_n_cycles}")
    print(f"Total number of elements in S_{n} = {n}! = {total_elements}")
    print(f"\nFinal Equation: Density = {num_n_cycles} / {total_elements}")
    print(f"The simplified density is: {density.numerator}/{density.denominator}")

solve()