from fractions import Fraction

# Number of i.i.d. random variables
n = 4

print("Step 1: Calculate E[Y]")
# Due to the symmetry of the problem, E[Y] = 1/2.
E_Y = Fraction(1, 2)
print(f"E[Y] = {E_Y}\n")

print("Step 2: Calculate E[Y^2]")
# Formula for the second moment of the k-th order statistic from n U[0,1] variables
# E[X_(k)^2] = k * (k + 1) / ((n + 1) * (n + 2))
def E_Xk_sq(k, n):
    return Fraction(k * (k + 1), (n + 1) * (n + 2))

# Calculate the second moments for k = 1, 2, 3, 4
E_X1_sq = E_Xk_sq(1, n)
E_X2_sq = E_Xk_sq(2, n)
E_X3_sq = E_Xk_sq(3, n)
E_X4_sq = E_Xk_sq(4, n)

print(f"The second moments of the order statistics for n={n} are:")
print(f"E[X_(1)^2] = {E_X1_sq}")
print(f"E[X_(2)^2] = {E_X2_sq}")
print(f"E[X_(3)^2] = {E_X3_sq}")
print(f"E[X_(4)^2] = {E_X4_sq}\n")

# Coefficients for the E[Y^2] formula based on the casework
c1 = Fraction(1, 8)
c2 = Fraction(3, 8)
c3 = Fraction(3, 8)
c4 = Fraction(1, 8)

# Calculate E[Y^2]
E_Y_sq = c1 * E_X1_sq + c2 * E_X2_sq + c3 * E_X3_sq + c4 * E_X4_sq

print("The formula for E[Y^2] is:")
print(f"E[Y^2] = {c1}*E[X_(1)^2] + {c2}*E[X_(2)^2] + {c3}*E[X_(3)^2] + {c4}*E[X_(4)^2]")
print("Plugging in the values:")
print(f"E[Y^2] = {c1}*{E_X1_sq} + {c2}*{E_X2_sq} + {c3}*{E_X3_sq} + {c4}*{E_X4_sq}")
print(f"E[Y^2] = {c1 * E_X1_sq} + {c2 * E_X2_sq} + {c3 * E_X3_sq} + {c4 * E_X4_sq}")
print(f"E[Y^2] = {E_Y_sq}\n")

print("Step 3: Calculate Var(Y)")
# Variance formula: Var(Y) = E[Y^2] - (E[Y])^2
Var_Y = E_Y_sq - E_Y**2

print("The formula for Var(Y) is: E[Y^2] - (E[Y])^2")
print(f"Var(Y) = {E_Y_sq} - ({E_Y})^2")
print(f"Var(Y) = {E_Y_sq} - {E_Y**2}")
print(f"Var(Y) = {Var_Y}\n")

print(f"The final calculated variance of Y is {Var_Y}, or approximately {float(Var_Y):.4f}.")
