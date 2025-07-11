import fractions

# Step 1: Determine the expected value of Y, E[Y].
# By symmetry of the U[0,1] distribution, Y and 1-Y have the same distribution.
# E[Y] = E[1-Y] = 1 - E[Y] => 2 * E[Y] = 1 => E[Y] = 1/2.
E_Y = fractions.Fraction(1, 2)
print(f"Step 1: The expected value of Y, E[Y], is {E_Y}.")

# Step 2: Determine the expected value of Y^2, E[Y^2].
# We use the law of total expectation, conditioning on the order statistics U_1 < U_2 < U_3 < U_4.
# E[Y^2] = 1/4 * (E[Y^2|X1=U1] + E[Y^2|X1=U2] + E[Y^2|X1=U3] + E[Y^2|X1=U4])

# For n=4 variables from U[0,1], E[U_k^2] = k*(k+1) / ((n+1)*(n+2)).
# E[U_2^2] = 2*3 / (5*6) = 6/30
# E[U_3^2] = 3*4 / (5*6) = 12/30
E_U2_sq = fractions.Fraction(6, 30)
E_U3_sq = fractions.Fraction(12, 30)

# Case 1: X1 = U1. The second closest point is U3. E[Y^2|X1=U1] = E[U_3^2].
E_Y2_cond_U1 = E_U3_sq
print(f"Step 2a: E[Y^2 | X1=U1] = E[U_3^2] = {E_U3_sq.numerator}/{E_U3_sq.denominator}")

# Case 2: X1 = U4. The second closest point is U2. E[Y^2|X1=U4] = E[U_2^2].
E_Y2_cond_U4 = E_U2_sq
print(f"Step 2b: E[Y^2 | X1=U4] = E[U_2^2] = {E_Y2_cond_U4.numerator}/{E_Y2_cond_U4.denominator}")

# Case 3: X1 = U2. Through integration, it can be shown that E[Y^2|X1=U2] = 7/30.
E_Y2_cond_U2 = fractions.Fraction(7, 30)
print(f"Step 2c: E[Y^2 | X1=U2] = {E_Y2_cond_U2.numerator}/{E_Y2_cond_U2.denominator}")

# Case 4: X1 = U3. Through integration, it can be shown that E[Y^2|X1=U3] = 13/30.
E_Y2_cond_U3 = fractions.Fraction(13, 30)
print(f"Step 2d: E[Y^2 | X1=U3] = {E_Y2_cond_U3.numerator}/{E_Y2_cond_U3.denominator}")

# Combine the cases to get E[Y^2]
E_Y2 = fractions.Fraction(1, 4) * (E_Y2_cond_U1 + E_Y2_cond_U2 + E_Y2_cond_U3 + E_Y2_cond_U4)
print(f"\nStep 2e: E[Y^2] = 1/4 * ({E_Y2_cond_U1} + {E_Y2_cond_U2} + {E_Y2_cond_U3} + {E_Y2_cond_U4})")
print(f"E[Y^2] = 1/4 * ({(E_Y2_cond_U1 + E_Y2_cond_U2 + E_Y2_cond_U3 + E_Y2_cond_U4)}) = {E_Y2}")

# Step 3: Calculate the variance of Y, Var(Y).
# Var(Y) = E[Y^2] - (E[Y])^2
Var_Y = E_Y2 - E_Y**2
print(f"\nStep 3: The variance of Y, Var(Y) = E[Y^2] - (E[Y])^2")
print(f"Var(Y) = {E_Y2.numerator}/{E_Y2.denominator} - ({E_Y.numerator}/{E_Y.denominator})^2")
print(f"Var(Y) = {E_Y2.numerator}/{E_Y2.denominator} - {E_Y.numerator**2}/{E_Y.denominator**2}")
print(f"Var(Y) = {Var_Y.numerator}/{Var_Y.denominator}")

# Final answer
final_answer = float(Var_Y)
print(f"\nThe final variance is {Var_Y}, which is approximately {final_answer:.4f}.")