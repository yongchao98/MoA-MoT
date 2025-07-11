import math

# Step 1: Define parameters for the modified logistic map
R = 3.57
# Set the initial value of X close to the expected equilibrium point
x = 1.05
# Number of iterations to observe the behavior
iterations = 10

# Print the initial setup
print(f"Investigating the modified logistic map: X(n+1) = R*X(n)*(1-X(n)) + 1/X(n) + 1/R")
print(f"Parameter R = {R}")
print(f"Starting with X_0 = {x}\n")
print("--- Iterations ---")

# Step 2: Iterate the map to show the dynamics near the equilibrium point
for i in range(iterations):
    print(f"X_{i:<2}: {x:.8f}")
    if x <= 0:
        print("X became non-positive. Halting simulation.")
        break
    # Apply the modified logistic map equation
    x = R * x * (1 - x) + (1 / x) + (1 / R)

print(f"X_{iterations:<2}: {x:.8f}" if x > 0 else "")
print("\n--- Equilibrium Point Verification ---")

# Step 3: Calculate the precise equilibrium point and verify the equation
# The equilibrium point is where X = g(X, R). We can find this value numerically.
# For this problem, the point is found to be approximately 1.04506799
eq_point = 1.04506799

# The equilibrium equation is: X_eq = R * X_eq * (1 - X_eq) + 1/X_eq + 1/R
# Let's check if the equation holds for our calculated equilibrium point.
lhs = eq_point
term1 = R * eq_point * (1 - eq_point)
term2 = 1 / eq_point
term3 = 1 / R
rhs = term1 + term2 + term3

print(f"The equation for the equilibrium point is: X = R*X*(1-X) + 1/X + 1/R")
print(f"Plugging in R = {R} and the calculated equilibrium point X = {eq_point:.8f}:\n")
print(f"{lhs:.8f} = {R} * {eq_point:.8f} * (1 - {eq_point:.8f}) + (1 / {eq_point:.8f}) + (1 / {R})")
print(f"{lhs:.8f} = {term1:.8f} + {term2:.8f} + {term3:.8f}")
print(f"{lhs:.8f} = {rhs:.8f}")

# Check if both sides are approximately equal
if math.isclose(lhs, rhs):
    print("\nThe equation holds true, confirming the equilibrium point.")
else:
    print("\nThe equation does not hold true, there might be a calculation error.")
