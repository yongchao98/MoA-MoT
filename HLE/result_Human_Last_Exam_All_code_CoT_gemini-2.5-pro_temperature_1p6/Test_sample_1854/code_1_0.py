import sympy

# Define the symbol for the constant c
c = sympy.Symbol('c')

# Based on the scaling analysis, the power of R on the left-hand side (LHS) of the inequality is derived.
# It comes from R^3 from scaling the integral and 1/R^2 from the norm of the operator T.
LHS_R_power = 3 - 2
print(f"The effective power of R on the LHS is {LHS_R_power}")

# The power of R on the right-hand side (RHS) comes from the problem statement R^{2c} (for the squared norm) and R^1 from scaling the norm.
# We ignore the epsilon for finding the boundary value of c.
RHS_R_power = 2*c + 1
print(f"The power of R on the RHS is {RHS_R_power}")

# To find the smallest c, we set the powers of R to be equal, representing the critical case.
equation = sympy.Eq(LHS_R_power, RHS_R_power)
print(f"\nSetting the powers equal gives the equation: {LHS_R_power} = {2*c + 1}")

# Solve the equation for c
solution = sympy.solve(equation, c)
c_value = solution[0]

print(f"The solution for c is: {c_value}")