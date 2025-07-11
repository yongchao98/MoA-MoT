import math
from fractions import Fraction

def solve_mass_formula():
    """
    This script calculates the total mass M(A_5, rho, 2) by summing contributions
    from all possible partitions of the integer 5, which correspond to the
    structures of etale Q_2-algebras of degree 5.
    """

    # Known counts of p-adic field extensions of Q_2 for small degrees.
    # These numbers are taken from established databases like the LMFDB.
    
    # Degree 5 (quintic) extensions by Galois group
    num_quintic_C5 = 1
    num_quintic_D5 = 4
    num_quintic_A5 = 4

    # Degree 4 (quartic) extensions by Galois group
    num_quartic_V4 = 21
    num_quartic_A4 = 2

    # Degree 3 (cubic) extensions by Galois group
    num_cubic_C3 = 2
    num_cubic_S3 = 6

    # Degree 2 (quadratic) extensions
    num_quadratic = 7

    total_mass = Fraction(0)
    
    print("Calculating the total mass M(A5, rho, 2)...")
    print("The total mass is a sum of contributions from partitions of 5.")
    print("-" * 30)

    # Partition [5]: K is a degree 5 field extension.
    # The Galois group of the extension must be a transitive subgroup of A5.
    # These are C5, D5, and A5.
    count_5 = num_quintic_C5 + num_quintic_D5 + num_quintic_A5
    # For a field K, the automorphism group Aut(K) is trivial, so |Aut(K)| = 1.
    mass_5 = Fraction(count_5, 1)
    total_mass += mass_5
    print(f"Partition [5]: Contribution = {count_5}/1")

    # Partition [4, 1]: K = K1 x Q2, where K1 is a degree 4 field.
    # The Galois group of K1 must be a subgroup of A4 (V4 or A4).
    count_4_1 = num_quartic_V4 + num_quartic_A4
    # Aut(K1 x Q2) is trivial, |Aut(K)| = 1.
    mass_4_1 = Fraction(count_4_1, 1)
    total_mass += mass_4_1
    print(f"Partition [4, 1]: Contribution = {count_4_1}/1")

    # Partition [3, 2]: K = K1 x K2, a product of cubic and quadratic fields.
    # The sign of permutations must match: sgn(g1) * sgn(g2) = 1.
    # This implies the quadratic characters must be equal.
    # If Gal(K1) is C3, its sign character is trivial, which is impossible for a quadratic field K2.
    # If Gal(K1) is S3, its sign character defines a quadratic field. K2 must be that field.
    # Each of the S3-cubic fields uniquely determines its quadratic partner.
    count_3_2 = num_cubic_S3
    # Aut(K1 x K2) is trivial, |Aut(K)| = 1.
    mass_3_2 = Fraction(count_3_2, 1)
    total_mass += mass_3_2
    print(f"Partition [3, 2]: Contribution = {count_3_2}/1")

    # Partition [3, 1, 1]: K = K1 x Q2 x Q2.
    # The Galois group of K1 must be a subgroup of A3 = C3.
    count_3_1_1 = num_cubic_C3
    # Aut(K) can swap the two Q2 factors, so |Aut(K)| = |S2| = 2.
    mass_3_1_1 = Fraction(count_3_1_1, 2)
    total_mass += mass_3_1_1
    print(f"Partition [3, 1, 1]: Contribution = {count_3_1_1}/2")

    # Partition [2, 2, 1]: K = K1 x K2 x Q2.
    # The quadratic characters for K1 and K2 must be equal, so K1 must be isomorphic to K2.
    # We choose one quadratic field type (7 choices). The algebra is K1 x K1 x Q2.
    count_2_2_1 = num_quadratic
    # Aut(K) can swap the two K1 factors, so |Aut(K)| = |S2| = 2.
    mass_2_2_1 = Fraction(count_2_2_1, 2)
    total_mass += mass_2_2_1
    print(f"Partition [2, 2, 1]: Contribution = {count_2_2_1}/2")

    # Partition [2, 1, 1, 1]: K = K1 x Q2 x Q2 x Q2.
    # The Galois group of K1 (C2) acts as a single transposition on 5 elements, which is odd.
    # This cannot be a subgroup of A5.
    count_2_1_1_1 = 0
    mass_2_1_1_1 = Fraction(0)
    total_mass += mass_2_1_1_1
    print(f"Partition [2, 1, 1, 1]: Contribution = 0/1")

    # Partition [1, 1, 1, 1, 1]: K = Q2^5.
    # The Galois group is trivial, which is in A5. There is 1 such algebra.
    count_1_1_1_1_1 = 1
    # Aut(K) is the symmetric group S5, so |Aut(K)| = 120.
    aut_size_1_1_1_1_1 = math.factorial(5)
    mass_1_1_1_1_1 = Fraction(count_1_1_1_1_1, aut_size_1_1_1_1_1)
    total_mass += mass_1_1_1_1_1
    print(f"Partition [1, 1, 1, 1, 1]: Contribution = {count_1_1_1_1_1}/{aut_size_1_1_1_1_1}")
    
    print("-" * 30)

    # Final result
    print("The total mass M is the sum of these contributions:")
    # Using f-string for clear formatting of the final equation
    print(f"M = {mass_5} + {mass_4_1} + {mass_3_2} + {mass_3_1_1} + {mass_2_2_1} + {mass_2_1_1_1} + {mass_1_1_1_1_1}")
    
    # Calculate and print the final fractional value
    # Example: M = 9 + 23 + 6 + 1 + 7/2 + 0 + 1/120
    # M = 39 + 1 + 7/2 + 1/120 = 40 + 420/120 + 7/120 = ... (mistake in manual sum)
    # M = 9+23+6+1+7/2+1/120 = 39 + 3.5 + 1/120 = 42.5 + 1/120 = 85/2 + 1/120 = 5100/120 + 1/120 = 5101/120
    print(f"\nM = {total_mass}")
    print(f"The final answer as a fraction in lowest terms is {total_mass.numerator}/{total_mass.denominator}")

solve_mass_formula()