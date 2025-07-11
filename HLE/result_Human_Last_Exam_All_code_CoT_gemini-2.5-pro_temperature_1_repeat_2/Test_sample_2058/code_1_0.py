import fractions

# This script calculates the total mass M(A_5, rho, 2)
# M = (1/|A_5|) * sum_{phi: Gamma_2 -> A_5} 2^(-c_rho(phi))
# The sum is decomposed by the image H of phi.
# M = (1/60) * sum_{H <= A_5} (# of H subgroups) * sum_{K: Gal(K/Q_2)=H} 2^(-c_rho(K))

# Contribution from H = C_1 (trivial group)
# 1 subgroup, 1 extension (Q_2), conductor c=0
s1 = 1 * fractions.Fraction(2, 1)**(-0)
print(f"Contribution from C1: {s1}")

# Contribution from H = C_2
# 15 subgroups. 7 extensions (1 unramified c=0, 2 tame c=2, 4 wild c=4)
s2 = 15 * (1 * fractions.Fraction(2, 1)**(-0) + 2 * fractions.Fraction(2, 1)**(-2) + 4 * fractions.Fraction(2, 1)**(-4))
print(f"Contribution from C2: {s2}")

# Contribution from H = C_3
# 10 subgroups. 2 extensions (1 unramified c=0, 1 tame c=2)
s3 = 10 * (1 * fractions.Fraction(2, 1)**(-0) + 1 * fractions.Fraction(2, 1)**(-2))
print(f"Contribution from C3: {s3}")

# Contribution from H = V_4
# 5 subgroups. Using the conductor formula a(K) which is the conductor exponent f from LMFDB
# 1 ext (f=0), 6 exts (f=2), 12 exts (f=3)
s4 = 5 * (1 * fractions.Fraction(2, 1)**(-0) + 6 * fractions.Fraction(2, 1)**(-2) + 12 * fractions.Fraction(2, 1)**(-3))
print(f"Contribution from V4: {s4}")

# Contribution from H = C_5
# 6 subgroups. 1 extension (unramified), c=0
s5 = 6 * (1 * fractions.Fraction(2, 1)**(-0))
print(f"Contribution from C5: {s5}")

# Contribution from H = D_3 (S_3)
# 10 subgroups. 7 extensions (1 unramified c=0, 6 ramified c=2)
s6 = 10 * (1 * fractions.Fraction(2, 1)**(-0) + 6 * fractions.Fraction(2, 1)**(-2))
print(f"Contribution from D3: {s6}")

# Contribution from H = D_5
# 6 subgroups. 1 extension (tame ramified), c=4
s7 = 6 * (fractions.Fraction(2, 1)**(-4))
print(f"Contribution from D5: {s7}")

# Contribution from H = A_4
# 5 subgroups. 4 extensions (1 unramified c=0, 3 ramified c=3)
s8 = 5 * (1 * fractions.Fraction(2, 1)**(-0) + 3 * fractions.Fraction(2, 1)**(-3))
print(f"Contribution from A4: {s8}")

# Total sum S
S = s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8
print(f"Total sum over homomorphisms is S = {s1} + {s2} + {s3} + {s4} + {s5} + {s6} + {s7} + {s8} = {S}")

# Total mass M = S / |A_5|
A5_order = 60
M = S / A5_order

# Final result as a fraction in lowest terms
final_fraction = fractions.Fraction(M).limit_denominator()
print(f"The total mass M(A_5, rho, 2) is S / |A_5| = {S} / {A5_order} = {M}")
print(f"The final answer as a fraction is: {final_fraction.numerator}/{final_fraction.denominator}")