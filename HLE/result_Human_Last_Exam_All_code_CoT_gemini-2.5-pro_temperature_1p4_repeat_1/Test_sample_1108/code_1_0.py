# The coefficients for the linearization of the predator-prey system
# at the non-trivial equilibrium point (1, 1).

# Coefficients of the Jacobian matrix A
a11 = -1
a12 = 1
a21 = -1
a22 = -1

# Coefficients of the constant vector B
b11 = 0
b22 = 0

# Print the values of the coefficients as requested
print("The coefficients of the linearized system are:")
print(f"a_11 = {a11}")
print(f"a_12 = {a12}")
print(f"a_21 = {a21}")
print(f"a_22 = {a22}")
print(f"b_11 = {b11}")
print(f"b_22 = {b22}")
