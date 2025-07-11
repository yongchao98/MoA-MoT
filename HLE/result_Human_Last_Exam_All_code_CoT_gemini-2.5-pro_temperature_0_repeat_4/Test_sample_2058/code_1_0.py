from fractions import Fraction

# This script calculates the total mass M(A_5, rho, 2).
# By a theorem of Kedlaya, this mass is equal to a_5(Q_2), the mass of quintic extensions of the 2-adic field Q_2.
# The mass a_n(Q_p) is the sum of 1/|Aut(K/Q_p)| over all degree n extensions K of Q_p.

print("Calculating the mass of quintic extensions of Q_2, denoted a_5(Q_2).")
print("Extensions are classified by ramification index 'e' and inertia degree 'f'. Here, ef=5.")
print("Since p=2 does not divide the degree 5, all extensions are tame.")
print("-" * 30)

# Case 1: Unramified extensions (e=1, f=5)
print("Case 1: Unramified extensions (e=1, f=5)")
# There is 1 unique unramified extension of degree 5.
num_unramified = 1
# It is a Galois extension, so its automorphism group has size 5.
aut_unramified = 5
mass_unramified = Fraction(num_unramified, aut_unramified)
print(f"Number of extensions: {num_unramified}")
print(f"Size of automorphism group: {aut_unramified}")
print(f"Mass contribution from this case: {mass_unramified.numerator}/{mass_unramified.denominator}")
print("-" * 30)

# Case 2: Totally ramified extensions (e=5, f=1)
print("Case 2: Totally ramified extensions (e=5, f=1)")
# The number of tame, totally ramified extensions of degree e is gcd(e, q-1), where q is the residue field size.
# For Q_2, q=2. So, number of extensions = gcd(5, 2-1) = 1.
num_ramified = 1
# This extension is not Galois, and since its degree is prime, its automorphism group is trivial.
aut_ramified = 1
mass_ramified = Fraction(num_ramified, aut_ramified)
print(f"Number of extensions: {num_ramified}")
print(f"Size of automorphism group: {aut_ramified}")
print(f"Mass contribution from this case: {mass_ramified.numerator}/{mass_ramified.denominator}")
print("-" * 30)

# Total Mass Calculation
total_mass = mass_unramified + mass_ramified
print("The total mass is the sum of the contributions from all cases.")
print(f"Total Mass = {mass_unramified.numerator}/{mass_unramified.denominator} + {mass_ramified.numerator}/{mass_ramified.denominator}")
print(f"Total Mass = {total_mass.numerator}/{total_mass.denominator}")

print(f"\nThus, the total mass M(A_5, rho, 2) is {total_mass.numerator}/{total_mass.denominator}.")