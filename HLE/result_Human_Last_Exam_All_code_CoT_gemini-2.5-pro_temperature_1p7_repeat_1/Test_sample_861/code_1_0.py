import numpy as np

# Step 1: Solve the characteristic equation 8r^2 - 6r - 1 = 0
# The roots are r = (6 ± sqrt(36 - 32)) / 16 = (6 ± 2) / 16
# So the roots are 8/16 = 1/2 and 4/16 = 1/4
r1 = 0.5
r2 = 0.25

# By convention, we assign the larger root to B and the smaller root to D.
B = r1
D = r2
print(f"The roots of the characteristic equation are {r1} and {r2}.")
print(f"Assigning the larger root to B and the smaller to D: B = {B}, D = {D}")

# Step 2: Find the particular solution constant E
# 8E - 6E + E = 1  => 3E = 1
E = 1/3
print(f"The particular solution is a constant E = {E:.4f}")

# Step 3: Use initial conditions to find A and C
# The general solution is y[n] = A*(B**n) + C*(D**n) + E
# We have a system of two linear equations from the initial conditions:
# y[0] = 1 => A*(B**0) + C*(D**0) + E = 1 => A + C = 1 - E
# y[-1] = 2 => A*(B**-1) + C*(D**-1) + E = 2 => A/B + C/D = 2 - E

# Let's set up the system in matrix form: M * x = v
# where x = [A, C]
M = np.array([
    [1, 1],
    [1/B, 1/D]
])

v = np.array([
    1 - E,
    2 - E
])

# Solve for A and C
try:
    coefficients = np.linalg.solve(M, v)
    A = coefficients[0]
    C = coefficients[1]
    print(f"Solving the system of equations yields: A = {A:.4f}, C = {C:.4f}")

    # Step 4: Print the full closed-form solution
    print("\nThe closed-form solution is:")
    print(f"y[n] = {A:.4f} * ({B})^n + {C:.4f} * ({D})^n + {E:.4f}")

    # Step 5: Calculate the final expression E/A + (D*C)/B
    result = E/A + (D*C)/B
    print("\nCalculating the expression E/A + (D*C)/B:")
    print(f"E/A = ({E:.4f}) / ({A:.4f}) = {E/A:.4f}")
    print(f"(D*C)/B = (({D}) * ({C:.4f})) / ({B}) = {(D*C)/B:.4f}")
    print(f"Final Result = {E/A:.4f} + {(D*C)/B:.4f} = {result:.4f}")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")

<<<0.75>>>