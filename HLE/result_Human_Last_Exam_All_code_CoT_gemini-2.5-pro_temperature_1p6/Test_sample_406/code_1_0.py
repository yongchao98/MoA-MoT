import numpy as np

# Define three functions F, G, H that satisfy the contraction conditions.
# For linear functions F(x,y,z) = c1*x + c2*y + c3*z + k, the condition is
# that the sum of the absolute values of the coefficients is less than 1.
# e.g., for F: |0.2| + |0.1| + |0.1| = 0.4 < 1.

def F(x, y, z):
    """Function F: X x Y x Z -> X"""
    # Contraction coefficients: 0.2, 0.1, 0.1. Constant: 5
    return 0.2 * x + 0.1 * y + 0.1 * z + 5

def G(y, z, x):
    """Function G: Y x Z x X -> Y"""
    # Contraction coefficients: 0.1, 0.2, 0.1. Constant: 2
    return 0.1 * y + 0.2 * z + 0.1 * x + 2

def H(z, x, y):
    """Function H: Z x X x Y -> Z"""
    # Contraction coefficients: 0.1, 0.1, 0.2. Constant: 3
    return 0.1 * z + 0.1 * x + 0.2 * y + 3

def find_tripled_fixed_point():
    """
    Iteratively finds the FGH-tripled fixed point.
    """
    # Initial guess for the fixed point (x, y, z)
    x, y, z = 0.0, 0.0, 0.0

    # The theorem guarantees convergence, so we iterate a fixed number of times.
    num_iterations = 100
    for _ in range(num_iterations):
        # Calculate the next point in the sequence
        x_next = F(x, y, z)
        y_next = G(y, z, x)
        z_next = H(z, x, y)
        
        # Update the current point
        x, y, z = x_next, y_next, z_next

    # The iteration converges to the fixed point (x, y, z).
    # Now we verify that F(x, y, z) = x, G(y, z, x) = y, and H(z, x, y) = z.
    
    # Calculate the results from the functions at the fixed point
    result_fx = F(x, y, z)
    result_gy = G(y, z, x)
    result_hz = H(z, x, y)
    
    # Print the verification equations. Due to floating point arithmetic,
    # the results will be extremely close but not perfectly identical.
    # We format the numbers to a reasonable precision.
    print("The FGH-tripled fixed point (x, y, z) is approximately ({:.4f}, {:.4f}, {:.4f})\n".format(x, y, z))
    print("Verifying the fixed point conditions:")
    
    # Print each number in the final equations
    print("F({:.4f}, {:.4f}, {:.4f}) = {:.4f} (approximates x)".format(x, y, z, result_fx))
    print("G({:.4f}, {:.4f}, {:.4f}) = {:.4f} (approximates y)".format(y, z, x, result_gy))
    print("H({:.4f}, {:.4f}, {:.4f}) = {:.4f} (approximates z)".format(z, x, y, result_hz))

if __name__ == "__main__":
    find_tripled_fixed_point()