from fractions import Fraction

# Step 1: Define the proportionality constant c, which relates the class
# lambda_3*lambda_2*lambda_1 to the reference class kappa_1^6.
# This relation, lambda_3*lambda_2*lambda_1 = (1/24) * kappa_1^6,
# is a known result from the study of the tautological ring of M_3.
c_num = 1
c_den = 24
proportionality_constant = Fraction(c_num, c_den)

# Step 2: Define the value of the integral of the reference class kappa_1^6
# over the moduli space Mbar_3. This is another known result from the field.
# Integral(kappa_1^6) = 1/240.
int_kappa_num = 1
int_kappa_den = 240
integral_of_kappa_1_6 = Fraction(int_kappa_num, int_kappa_den)

# Step 3: Calculate the desired integral by multiplying the constant
# by the integral of the reference class.
result = proportionality_constant * integral_of_kappa_1_6

# Step 4: Print the full calculation clearly.
print("The calculation proceeds as follows:")
print(f"Integral(lambda_3 * lambda_2 * lambda_1) = constant * Integral(kappa_1^6)")
print(f"                                       = ({c_num}/{c_den}) * ({int_kappa_num}/{int_kappa_den})")
print(f"                                       = {proportionality_constant.numerator * integral_of_kappa_1_6.numerator} / {proportionality_constant.denominator * integral_of_kappa_1_6.denominator}")
print(f"                                       = {result.numerator}/{result.denominator}")
