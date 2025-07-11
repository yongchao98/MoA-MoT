import fractions

# M(A_5, 1, 2)
# Using the known value |Hom(Gamma_2, A_5)| = 742 and |A_5| = 60
m1_num = 742
m1_den = 60
m1 = fractions.Fraction(m1_num, m1_den)
print(f"M(A_5, 1, 2) = {m1_num}/{m1_den} = {m1}")

# M(A_5, chi_4, 2)
# Sum over homomorphisms = unramified_contrib + ramified_contrib
# Unramified contribution to the sum
unram_contrib = 60

# Ramified contribution from C2 images
c2_contrib_f2 = fractions.Fraction(30, 16) # 2 fields * 15 subgroups * 2^(-2*2)
c2_contrib_f3 = fractions.Fraction(60, 64) # 4 fields * 15 subgroups * 2^(-2*3)
c2_contrib = c2_contrib_f2 + c2_contrib_f3

# Ramified contribution from V4 images
v4_contrib_d6 = fractions.Fraction(90, 64) # 3 fields * 30 maps * 2^(-6)
v4_contrib_d8 = fractions.Fraction(120, 256) # 4 fields * 30 maps * 2^(-8)
v4_contrib = v4_contrib_d6 + v4_contrib_d8

# Ramified contribution from S3 images
s3_contrib = fractions.Fraction(60, 4) # 60 maps * 2^(-2)

# Total sum for chi_4
total_sum_chi4 = unram_contrib + c2_contrib + v4_contrib + s3_contrib
print(f"Sum for chi_4 = {unram_contrib} + {c2_contrib} + {v4_contrib} + {s3_contrib} = {total_sum_chi4}")


# M(A_5, chi_4, 2) = (1/|A_5|) * Sum
m4_den = 60
m4 = total_sum_chi4 / m4_den
print(f"M(A_5, chi_4, 2) = {total_sum_chi4} / {m4_den} = {m4}")


# Total Mass M(A_5, rho, 2) = M(A_5, 1, 2) + M(A_5, chi_4, 2)
total_mass = m1 + m4
print(f"\nTotal Mass M(A_5, rho, 2) = M(A_5, 1, 2) + M(A_5, chi_4, 2)")
print(f"M(A_5, rho, 2) = {m1} + {m4}")
print(f"M(A_5, rho, 2) = {m1.numerator*m4.denominator + m4.numerator*m1.denominator}/{m1.denominator*m4.denominator} = {total_mass}")

print(f"\nThe equation is {m1.numerator}/{m1.denominator} + {m4.numerator}/{m4.denominator} = {total_mass.numerator}/{total_mass.denominator}")
