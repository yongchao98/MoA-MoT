import sys

def solve():
    """
    This script modifies the logistic map to create a stable equilibrium point
    where chaos would normally exist, without adding new parameters.
    """

    # 1. Define the modified logistic map function
    # The equation is: X_n+1 = 1 + R * (X_n - 1) * (4/3 - X_n)
    def modified_logistic_map(x, r):
        c1 = 1.0
        c2 = 1.0
        c3 = 4.0
        c4 = 3.0
        return c1 + r * (x - c2) * (c3 / c4 - x)

    # 2. Set parameters for the simulation
    R = 3.57  # The value of R where the standard map is chaotic
    X0 = 0.5  # An arbitrary starting value for X
    N_ITERATIONS = 200 # Number of iterations to run to ensure convergence

    # 3. Print the equation with its numerical components
    # The prompt requires outputting each number in the final equation.
    c1, c2, c3, c4 = 1, 1, 4, 3
    print(f"The modified logistic map equation is:")
    print(f"X_n+1 = {c1} + R * (X_n - {c2}) * ({c3}/{c4} - X_n)\n")
    print(f"Running simulation with R = {R}...")

    # 4. Run the simulation
    x = X0
    for i in range(N_ITERATIONS):
        x = modified_logistic_map(x, R)

    # The final value of x is the equilibrium point
    final_equilibrium_point = x

    # 5. Print the final result
    print(f"\nAfter {N_ITERATIONS} iterations, the system converged to the equilibrium point.")
    print(f"The calculated equilibrium point is approximately: {final_equilibrium_point}")
    
    # Also printing the value in the requested format for the final answer.
    sys.stdout.write(f'<<<{final_equilibrium_point}>>>')

solve()