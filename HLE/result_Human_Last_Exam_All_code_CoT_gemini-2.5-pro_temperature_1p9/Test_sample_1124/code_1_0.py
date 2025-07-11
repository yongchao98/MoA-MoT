import math

# Plan:
# 1. Define a modified logistic map to create a stable equilibrium point at X=1.05 when R=3.57.
#    The standard logistic map X_n+1 = R*X_n*(1-X_n) is chaotic at R=3.57.
# 2. To ensure stability, we need the derivative of our new map g(X) to be less than 1 at the equilibrium point.
#    A superstable point (derivative is 0) can be created with a map of the form: g(X) = X* + C*(X-X*)^2.
# 3. We set the equilibrium point X* to 1.05. We can set C = R to incorporate R into the map.
#    This gives the final equation: X_n+1 = 1.05 + R * (X_n - 1.05)^2.
# 4. This script will demonstrate that this new map converges to 1.05.

# --- Parameters ---
R = 3.57
X_STAR = 1.05
INITIAL_X = 0.5
ITERATIONS = 20

def modified_logistic_map(x, r):
    """
    This is the modified logistic map designed to have a superstable
    fixed point at X_STAR.
    """
    return X_STAR + r * (x - X_STAR)**2

# --- Main Execution ---
print("This script demonstrates a modified logistic map.\n")

print("The final equation for the modified map is of the form:")
print("X_n+1 = C1 + R * (X_n - C2)^2\n")

# Output each number in the final equation
print(f"For this problem, the numbers in the equation are:")
print(f"C1 = {X_STAR}")
print(f"C2 = {X_STAR}")
print(f"R  = {R} (the value where the standard map is chaotic)\n")

print(f"The specific equation being simulated is:")
print(f"X_n+1 = {X_STAR} + {R} * (X_n - {X_STAR})^2\n")


# Simulate the map to show convergence
x = INITIAL_X
print(f"Simulating for {ITERATIONS} iterations, starting with X_0 = {x}...\n")
# Print header for the results table
print(f"{'Iteration':<10} | {'Value of X':<20}")
print("-" * 33)
print(f"{0:<10} | {x:<20.15f}")

for i in range(ITERATIONS):
    x = modified_logistic_map(x, R)
    print(f"{i+1:<10} | {x:<20.15f}")
    # Stop if it has converged to avoid spamming the same number
    if math.isclose(x, X_STAR):
        print("\nValue has converged to the stable equilibrium point.")
        break