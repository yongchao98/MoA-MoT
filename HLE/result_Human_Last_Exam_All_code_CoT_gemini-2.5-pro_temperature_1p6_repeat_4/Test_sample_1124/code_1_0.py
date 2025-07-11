import sys

def modified_logistic_map(x, R):
    """
    A modified logistic map based on the Beverton-Holt model.
    The constant C is derived from the problem's conditions.
    """
    # The constant C is derived to set the equilibrium point at 1.05 for R=3.57
    # C = R_chaos - X_target = 3.57 - 1.05 = 2.52
    C = 2.52
    return R * x / (C + x)

# Set the parameters from the problem
R_chaos = 3.57
target_equilibrium = 1.05
C_derived = R_chaos - target_equilibrium

# Initial condition for the iteration
x0 = 0.5 
n_iterations = 15

# --- Output the results ---

# Print the final equation as requested
print("The modified map equation is:")
print("X_n+1 = {} * X_n / ({} + X_n)".format(R_chaos, round(C_derived, 2)))
print("-" * 40)

# Iterate the map and show convergence
x = x0
print("Iterating the map starting from X_0 = {}".format(x0))
print("X_0: {:.6f}".format(x))
for i in range(n_iterations):
    x = modified_logistic_map(x, R_chaos)
    print("X_{}: {:.6f}".format(i + 1, x))

print("\nThe system converges to the stable equilibrium point.")

# Final answer format
print("Final equilibrium point is approximately: {}".format(round(x, 2)))
sys.stdout.write("<<<1.05>>>")