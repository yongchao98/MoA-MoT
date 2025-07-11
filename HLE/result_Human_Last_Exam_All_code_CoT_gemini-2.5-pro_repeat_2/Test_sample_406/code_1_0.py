import math

# Define three functions F, G, H that satisfy the contraction condition.
# For these linear functions, the sum of the absolute values of the coefficients
# is 0.1+0.1+0.1 + 0.2+0.1+0.05 + 0.1+0.1+0.05 = 0.9, which is less than 1.
# This guarantees that a unique tripled fixed point exists.
def F(x, y, z):
    """Function F: R*R*R -> R"""
    return 0.1 * x + 0.1 * y - 0.1 * z + 5

def G(y, z, x):
    """Function G: R*R*R -> R"""
    return 0.2 * y + 0.1 * z + 0.05 * x + 2

def H(z, x, y):
    """Function H: R*R*R -> R"""
    return -0.1 * z + 0.1 * x - 0.05 * y + 1

def find_tripled_fixed_point():
    """
    Iteratively finds the tripled fixed point (x, y, z) such that:
    F(x, y, z) = x
    G(y, z, x) = y
    H(z, x, y) = z
    """
    # Initial guess
    x, y, z = 0.0, 0.0, 0.0

    # Iteration parameters
    max_iterations = 1000
    tolerance = 1e-9

    print("Searching for the FGH-tripled fixed point...")
    print("Initial point (x,y,z) = (0.0, 0.0, 0.0)\n")

    for i in range(max_iterations):
        # Calculate the next point in the sequence
        x_new = F(x, y, z)
        y_new = G(y, z, x)
        z_new = H(z, x, y)

        # Check for convergence
        if (math.isclose(x, x_new, rel_tol=tolerance) and
            math.isclose(y, y_new, rel_tol=tolerance) and
            math.isclose(z, z_new, rel_tol=tolerance)):
            print(f"Convergence reached after {i+1} iterations.\n")
            x, y, z = x_new, y_new, z_new
            break

        # Update the current point
        x, y, z = x_new, y_new, z_new
    else:
        print("Did not converge within the maximum number of iterations.\n")

    # Print the results in the required format
    print(f"The calculated tripled fixed point is (x, y, z) = ({x:.4f}, {y:.4f}, {z:.4f})\n")
    print("--- Final Equations ---")
    
    # Output each number in the final equation
    print(f"F({x:.4f}, {y:.4f}, {z:.4f}) = {F(x, y, z):.4f}, which is approximately x ({x:.4f})")
    print(f"G({y:.4f}, {z:.4f}, {x:.4f}) = {G(y, z, x):.4f}, which is approximately y ({y:.4f})")
    print(f"H({z:.4f}, {x:.4f}, {y:.4f}) = {H(z, x, y):.4f}, which is approximately z ({z:.4f})")


if __name__ == '__main__':
    find_tripled_fixed_point()
