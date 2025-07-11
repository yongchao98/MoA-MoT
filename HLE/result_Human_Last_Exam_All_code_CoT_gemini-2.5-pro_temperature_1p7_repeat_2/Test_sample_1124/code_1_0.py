import math

def solve():
    """
    This function implements the modified logistic map, demonstrates its convergence
    to a stable equilibrium point, and prints the results as per the problem description.
    """
    # The value of R where the standard logistic map is chaotic
    R = 3.57
    
    # An arbitrary starting point for X
    x = 0.5

    # Define the function a(R)
    a_R = math.e - 1/R
    
    # Calculate the theoretical equilibrium point X*
    # X* = (R - 1) / a(R)
    x_equilibrium = (R - 1) / a_R

    # --- Output Section ---
    print("Modified Map Equation: X_n+1 = R * X_n / (1 + a(R) * X_n)")
    print("where a(R) = e - 1/R")
    print("-" * 40)
    print(f"For R = {R}:")
    print(f"a({R}) = e - 1/{R} = {math.e:.5f} - {1/R:.5f} = {a_R:.5f}")
    print(f"Theoretical equilibrium point X* = ({R} - 1) / {a_R:.5f} = {x_equilibrium:.5f}")
    print("-" * 40)

    print(f"Simulating convergence from X_0 = {x}:")
    # Simulate for a number of iterations to show convergence
    n_iterations = 50
    for i in range(n_iterations):
        x_next = R * x / (1 + a_R * x)
        
        # Print the first few and the final iteration to show the numbers
        if i < 5 or i == n_iterations - 1:
            # Outputting each number in the final equation for the current step
            print(f"Step {i+1}: X_{i+1} = {R} * {x:.5f} / (1 + {a_R:.5f} * {x:.5f}) = {x_next:.5f}")
        x = x_next

    print("-" * 40)
    print(f"The simulation converges to a stable point approximately equal to {x:.5f}")

solve()
<<<1.05407>>>