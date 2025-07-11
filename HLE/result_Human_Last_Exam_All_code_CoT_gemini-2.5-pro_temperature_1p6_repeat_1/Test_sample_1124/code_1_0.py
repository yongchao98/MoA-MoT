import math

def run_modified_logistic_map():
    """
    This function implements a modified logistic map that creates a stable
    equilibrium point instead of chaos for R = 3.57.
    """
    # Define the R value specified in the problem
    R = 3.57

    # Set an initial value for X
    x = 0.5

    # Set the number of iterations to demonstrate convergence
    iterations = 100

    # This is the modified logistic map function based on the Ricker model.
    # The growth factor K is defined in terms of R as K = R - 1 + 1/R.
    def modified_map(x_n, r_val):
        k_term = r_val - 1 + (1 / r_val)
        return k_term * x_n * math.exp(-x_n)

    print(f"Simulating the modified logistic map for R = {R}")
    print("This map should converge to a stable equilibrium point.")
    print("-" * 50)

    # Iterate the map to find the equilibrium point
    for i in range(iterations):
        x = modified_map(x, R)
        # Print the first few steps to show the convergence process
        if i < 15:
            print(f"Iteration {i+1:2d}: X = {x:.6f}")
        elif i == 15:
            print("...")

    print(f"\nAfter {iterations} iterations, the system has stabilized.")
    final_x = x
    print(f"The final equilibrium point is: {final_x:.6f}")

    # The problem asks to output the numbers in the final equation.
    # At equilibrium, X_n+1 = X_n. We show the calculation for the final step.
    k_val = R - 1 + (1 / R)
    print("\nThe final state equation with the numbers filled in:")
    print(f"X_n+1 = (R - 1 + 1/R) * X_n * exp(-X_n)")
    print(f"{modified_map(final_x, R):.6f} = ({R} - 1 + {1/R:.6f}) * {final_x:.6f} * exp(-{final_x:.6f})")

    # For verification, we can calculate the theoretical fixed point: ln(K)
    theoretical_x = math.log(k_val)
    print(f"\nThe theoretical equilibrium point is ln({k_val:.6f}) = {theoretical_x:.6f}")


if __name__ == "__main__":
    run_modified_logistic_map()
