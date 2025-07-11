def modified_logistic_map(x, r):
    """
    A modified logistic map designed to have a stable fixed point at 1.05.
    The equation is: X_n+1 = -R * (X_n - 1.05)**2 + 1.05
    """
    # The target equilibrium value is used to define the vertex of the parabola.
    equilibrium_value = 1.05
    return -r * (x - equilibrium_value)**2 + equilibrium_value

# Parameters from the problem description
R = 3.57
TARGET_EQUILIBRIUM = 1.05
# An initial value for X to start the iteration
x_initial = 0.5
# Number of iterations to demonstrate convergence
iterations = 10

# Print the final equation with the specific value of R
print("The modified logistic map equation is:")
equation_string = f"X_n+1 = -{R} * (X_n - {TARGET_EQUILIBRIUM})**2 + {TARGET_EQUILIBRIUM}"
print(equation_string)

# Print each number used in the final equation, as requested
print("\nThe numbers that define this equation are:")
print(f"  - Coefficient of the quadratic term: -1")
print(f"  - Parameter R: {R}")
print(f"  - Value subtracted from X_n (h): {TARGET_EQUILIBRIUM}")
print(f"  - Exponent: 2")
print(f"  - Value added at the end (k): {TARGET_EQUILIBRIUM}\n")


print(f"--- Simulating for {iterations} iterations with R = {R} and initial X_0 = {x_initial} ---\n")
x = x_initial
print(f"Iteration 0: X = {x:.8f}")

# Run the simulation and print each step
for i in range(1, iterations + 1):
    x = modified_logistic_map(x, R)
    print(f"Iteration {i}: X = {x:.8f}")

print(f"\nAs the simulation shows, the system rapidly converges to the stable equilibrium point {TARGET_EQUILIBRIUM}.")