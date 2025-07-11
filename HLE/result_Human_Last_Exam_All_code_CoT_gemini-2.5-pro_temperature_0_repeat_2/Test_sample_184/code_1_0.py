import math

def solve_sum():
    """
    Calculates the sum of 1/n^2 for all positive integers n where no prime
    divides n with multiplicity 1, 2, or 5.
    """

    # The problem asks to evaluate the sum S = Σ_{n∈S} 1/n^2, where S is the set of
    # positive integers n such that for any prime p dividing n, its exponent is not 1, 2, or 5.
    # This sum can be expressed as an Euler product over all primes p:
    # S = Π_p [ Σ_{k∈A} (1/p²)^k ], where A = {0, 3, 4, 6, 7, ...} is the set of allowed exponents.

    # The key insight is that the set of allowed exponents A is the set of all numbers
    # that can be written as 3a + 4b for non-negative integers a and b.
    # The numbers that cannot be formed this way are 1, 2, and 5, which are exactly the forbidden exponents.

    # This allows us to rewrite the inner sum for each prime p. Let x = 1/p²:
    # Σ_{k∈A} x^k = Σ_{a≥0, b≥0} x^(3a+4b) = (Σ_{a≥0} (x³)^a) * (Σ_{b≥0} (x⁴)^b)
    # This simplifies to (1/(1-x³)) * (1/(1-x⁴)).

    # Substituting x = 1/p² back, the term for each prime is (1/(1-1/p⁶)) * (1/(1-1/p⁸)).
    # The total sum is the product over all primes:
    # S = [ Π_p (1/(1-1/p⁶)) ] * [ Π_p (1/(1-1/p⁸)) ]

    # These two products are the definitions of the Riemann zeta function ζ(6) and ζ(8).
    # S = ζ(6) * ζ(8)

    # The known values for these zeta functions are:
    # ζ(6) = π⁶ / 945
    # ζ(8) = π⁸ / 9450

    zeta6_den = 945
    zeta8_den = 9450

    final_den = zeta6_den * zeta8_den
    final_pi_power = 6 + 8

    print("The sum can be expressed as the product of two Riemann zeta functions: ζ(6) * ζ(8).")
    print(f"The value of ζ(6) is π^6 / {zeta6_den}.")
    print(f"The value of ζ(8) is π^8 / {zeta8_den}.")
    print("\nThe final sum is the product of these two values:")
    print(f"Sum = (π^6 / {zeta6_den}) * (π^8 / {zeta8_den})")
    print(f"    = π^({6 + 8}) / ({zeta6_den} * {zeta8_den})")
    print(f"    = π^{final_pi_power} / {final_den}")

solve_sum()
<<<π^14 / 8930250>>>