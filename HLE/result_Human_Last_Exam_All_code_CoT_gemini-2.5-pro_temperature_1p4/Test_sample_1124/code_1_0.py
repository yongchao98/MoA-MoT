import math

def modified_logistic_map(r, x, iterations=100):
    """
    Iterates the modified logistic map to find the equilibrium point.
    The modified map is: X_{n+1} = X_n * (1 - X_n) + R / 4
    """
    for _ in range(iterations):
        # This check is to prevent overflow if x becomes too large or small
        if not -1e6 < x < 1e6:
            return float('inf') 
        x = x * (1 - x) + r / 4
    return x

# System parameters
R = 3.57
X0 = 0.5  # Initial value for X

# Calculate the equilibrium point by iterating the map
equilibrium_point = modified_logistic_map(R, X0)

# The equilibrium point is derived from the equation: X^2 = R / 4
# Which simplifies to X = sqrt(R) / 2
print("The modified logistic map is: X_n+1 = X_n * (1 - X_n) + R / 4")
print(f"For R = {R}, the system converges to a stable equilibrium point.")
print("\nThe equilibrium equation is: X = sqrt(R) / 2")

# Print the final equation with the numbers plugged in
# The format is: result = equation_with_numbers
print(f"{equilibrium_point:.4f} = sqrt({R}) / 2")