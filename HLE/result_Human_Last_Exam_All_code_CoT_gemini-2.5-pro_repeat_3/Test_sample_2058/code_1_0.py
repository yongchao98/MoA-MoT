from fractions import Fraction

# The contributions from each solvable subgroup H of A_5 to the sum part of the mass formula.
# S_H = sum over extensions K with Galois group H of |Inj(H, A_5)| * 2^(-c(psi_4|_H, K))
S_C1 = Fraction(1)
S_C2 = Fraction(285, 16)
S_C3 = Fraction(25)
S_V4 = Fraction(105, 32)
S_C5 = Fraction(0)
S_S3 = Fraction(555, 32)
S_D5 = Fraction(60)
S_A4 = Fraction(15, 4)

# The total sum over all homomorphisms
total_sum = S_C1 + S_C2 + S_C3 + S_V4 + S_C5 + S_S3 + S_D5 + S_A4

# The order of the group A_5
order_A5 = 60

# The total mass M(A_5, rho, 2)
M = Fraction(1, order_A5) * total_sum

# The equation for the total mass is:
# M = (1/|A_5|) * (S_C1 + S_C2 + S_C3 + S_V4 + S_C5 + S_S3 + S_D5 + S_A4)
print("The calculation for the total mass M(A_5, rho, 2) is:")
print(f"M = (1/{order_A5}) * ({S_C1} + {S_C2} + {S_C3} + {S_V4} + {S_C5} + {S_S3} + {S_D5} + {S_A4})")
print(f"M = (1/{order_A5}) * ({total_sum})")
print(f"M = {M}")
print(f"The total mass is {M.numerator}/{M.denominator}.")
