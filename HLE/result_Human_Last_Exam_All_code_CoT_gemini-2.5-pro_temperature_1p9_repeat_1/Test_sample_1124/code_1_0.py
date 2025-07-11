def modified_logistic_map(x, r):
    """
    This is the modified logistic map function.
    X_n+1 = R * X_n * (1 - X_n) + X_n^R
    """
    term1 = r * x
    term2 = (1 - x)
    term3 = x**r
    x_next = term1 * term2 + term3
    return x_next, term1, term2, term3

# Define the parameters for the map
R = 3.57
# Set the initial value for X
x_n = 0.5
# Set the number of iterations to observe
iterations = 10

# Explain the modified equation
print("The modified logistic map is: X_n+1 = R * X_n * (1 - X_n) + X_n^R")
print(f"\nRunning with R = {R} and starting with X_0 = {x_n}\n")

# Iterate through the map to show convergence
for i in range(iterations):
    # Get the next value and the components of the calculation
    x_n_plus_1, term1_val, term2_val, term3_val = modified_logistic_map(x_n, R)

    # Print the equation with the specific numbers for the current step
    # This shows each number used in the final equation as requested
    print(f"Step {i}:")
    print(f"  X_{i+1} = {R} * {x_n:.6f} * (1 - {x_n:.6f}) + ({x_n:.6f})^{R}")
    
    # Break down the calculation
    calc_term1 = term1_val * term2_val
    print(f"        = {calc_term1:.6f} + {term3_val:.6f}")
    
    # Show the final result for the step
    print(f"        = {x_n_plus_1:.6f}\n")
    
    # Update x_n for the next iteration
    x_n = x_n_plus_1
    
    # A check to stop early if convergence is reached
    if abs(x_n - 1.0) < 1e-9:
        print("The system has converged to the equilibrium point.")
        break

print(f"The final equilibrium point is: {x_n:.6f}")