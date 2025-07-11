# The modified logistic map equation is X_n+1 = R * X_n * (1 - X_n) + X_n^3

# Parameters for the simulation
R = 3.57
# Initial value for X, starting in the typical range for the standard map
X = 0.5
# Number of iterations to run the map to ensure convergence
iterations = 100

# The numbers in the final equation are R=3.57, the constant 1, and the exponent 3.
# The equation is constructed using these numbers.
print("The modified logistic map equation is:")
print(f"X_n+1 = {R} * X_n * (1 - X_n) + X_n^3")
print("-" * 30)

# Iterate the map to find the equilibrium point
for i in range(iterations):
    X = R * X * (1 - X) + X**3

print(f"For R = {R}, the system was iterated {iterations} times.")
print("Instead of chaos, it converges to a stable equilibrium point at approximately:")
# The result should be very close to 1.0
print(X)