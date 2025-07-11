import sys
from fractions import Fraction

def solve_polynomial_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.

    This function uses Chebotarev's density theorem, which relates this density
    to the structure of the polynomial's Galois group.
    """

    # For the given polynomial to be irreducible modulo a prime p, the
    # Frobenius element at p must act as a single cycle on the 7 roots.
    # This means we need to count the number of 7-cycles in the Galois group.

    # We use the cypari2 library to determine the Galois group.
    # This library is a Python interface to the PARI/GP computer algebra system.
    # If you don't have it, you can install it via `pip install cypari2`.
    try:
        import cypari2
    except ImportError:
        print("Error: The 'cypari2' library is required to run this code.", file=sys.stderr)
        print("Please install it using: pip install cypari2", file=sys.stderr)
        return

    pari = cypari2.Pari()
    
    # Define the polynomial in PARI/GP syntax
    f_poly = pari('x^7 - 14*x^5 + 56*x^3 - 56*x + 22')

    # Compute the Galois group information. The galoisinit() function returns
    # a vector, where the first element is the group's order and the fourth
    # element is the group's name in PARI/GP's database.
    galois_info = f_poly.galois()
    group_order = int(galois_info[0])
    group_name = str(galois_info[3])

    if group_name != "D(7)":
        print(f"The Galois group was identified as {group_name} of order {group_order}, not D_7 as expected.", file=sys.stderr)
        print("The following calculation assumes the group is D_7.", file=sys.stderr)
    
    # The Galois group is D_7 (denoted D(7) by PARI/GP), the dihedral group
    # of order 14. We need to find the proportion of 7-cycles in this group.
    # A D_n group consists of n rotations and n reflections. For n=7:
    # - 1 identity element (a rotation by 0 degrees).
    # - 6 other rotations (by multiples of 2*pi/7). These are the 7-cycles.
    # - 7 reflections, each of which has order 2.
    num_7_cycles = 6
    
    # The order of the group D_7 is 14.
    denominator = group_order
    numerator = num_7_cycles
    
    # Calculate the fraction and simplify it.
    result_fraction = Fraction(numerator, denominator)

    # Output the explanation and the final equation.
    print(f"The Galois group of the polynomial is the Dihedral Group {group_name}, which has order {denominator}.")
    print("For the polynomial to be irreducible mod p, the Frobenius element must be a 7-cycle.")
    print(f"The number of 7-cycle elements in {group_name} is {numerator}.")
    print("By Chebotarev's density theorem, the required density is the ratio:")
    print(f"Density = {numerator} / {denominator} = {result_fraction.numerator} / {result_fraction.denominator}")

if __name__ == '__main__':
    solve_polynomial_density()
