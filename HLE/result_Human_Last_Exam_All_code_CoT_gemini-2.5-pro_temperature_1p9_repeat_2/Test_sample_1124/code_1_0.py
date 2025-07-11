import sys

# Step 1: Define the given parameter R and the target equilibrium.
# The standard logistic map X_n+1 = R * X_n * (1 - X_n) is chaotic at this R.
R = 3.57
# The target equilibrium point is approximately 1.05.

# Step 2: Formulate the modified logistic map.
# We found a relationship R / 1.05 = 3.4.
# We will construct a new map where the equilibrium X_eq is exactly R / 3.4.
# The modified map is: X_n+1 = X_n * (1 + R/C - X_n) where C is our constant.
C = 3.4

# The equilibrium point for this map is found by setting X_n+1 = X_n = X_eq:
# X_eq = X_eq * (1 + R/C - X_eq)
# Dividing by X_eq (for the non-trivial solution) gives:
# 1 = 1 + R/C - X_eq
# Which solves to: X_eq = R/C

# Calculate the equilibrium point for the given R.
eq_point = R / C

# The derivative of our map f(X) = X * (1 + R/C - X) is f'(X) = 1 + R/C - 2*X.
# At the equilibrium point X_eq = R/C, the derivative is:
# f'(X_eq) = 1 + R/C - 2*(R/C) = 1 - R/C.
# The stability condition is |f'(X_eq)| < 1.
# For R=3.57, f'(X_eq) = 1 - 3.57/3.4 = 1 - 1.05 = -0.05.
# Since |-0.05| < 1, the equilibrium point is stable.

# Step 3: Print the final equation and demonstrate its behavior.
print("--- Modified Logistic Map ---")
print("The standard logistic map is X_n+1 = R * X_n * (1 - X_n)")
print(f"For R = {R}, we want an equilibrium point of approximately 1.05.")
print("\nThe modified equation that satisfies this is:")
# Flush the buffer to ensure the explanation is printed before any potential errors
sys.stdout.flush()

print("\n\tX_n+1 = X_n * (1 + R/C - X_n)\n")
print("Where the numbers in the final equation are:")
print(f"\tConstant '1'")
print(f"\tConstant 'C' = {C}")
print(f"\tParameter 'R' = {R}")

print("\nThis gives a stable equilibrium point at X_eq = R / C")
print(f"X_eq = {R} / {C} = {eq_point}")

print("\n--- Simulation ---")
print("Starting with X_0 = 0.5, the system will converge to the new equilibrium:")
# Start with an initial value
x = 0.5
print(f"X_0 = {x}")

# Iterate a few times to show convergence
for i in range(1, 10):
    x = x * (1 + R/C - x)
    print(f"X_{i} = {x:.6f}")
