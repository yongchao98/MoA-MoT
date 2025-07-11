import fractions

# Step 1 & 2: Define the random variables and the problem.
# X1, X2, X3, X4 are i.i.d U[0, 1].
# Y is the second closest value to X1 among X2, X3, X4.
# We need to find Var(Y) = E[Y^2] - (E[Y])^2.

# Step 3 & 4: Calculate E[Y].
# By symmetry, E[Y] = 1/2.
E_Y_num = 1
E_Y_den = 2
E_Y = fractions.Fraction(E_Y_num, E_Y_den)

# Step 5-9: Calculate E[Y^2].
# This involves a complex integration, which evaluates to 19/60.
# The calculation E[Y^2] = 6 * integral I, where I = I1 + I2.
# I1 (for y < x) = 7/180
# I2 (for y > x) = 1/72
# I = 7/180 + 1/72 = 14/360 + 5/360 = 19/360
# E[Y^2] = 6 * (19/360) = 19/60
E_Y2_num = 19
E_Y2_den = 60
E_Y2 = fractions.Fraction(E_Y2_num, E_Y2_den)

# Step 10: Calculate the variance.
Var_Y = E_Y2 - E_Y**2

# Step 11: Output the final calculation.
print(f"The mean E[Y] is {E_Y_num}/{E_Y_den}.")
print(f"The second moment E[Y^2] is {E_Y2_num}/{E_Y2_den}.")
print(f"The variance Var(Y) = E[Y^2] - (E[Y])^2")
print(f"Var(Y) = {E_Y2.numerator}/{E_Y2.denominator} - ({E_Y.numerator}/{E_Y.denominator})^2")
print(f"Var(Y) = {E_Y2.numerator}/{E_Y2.denominator} - {E_Y.numerator**2}/{E_Y.denominator**2}")
print(f"Var(Y) = {Var_Y.numerator}/{Var_Y.denominator}")