import math

def F(x, y, z):
    """
    Sample function F: R*R*R -> R.
    This function is part of a system that satisfies contraction conditions.
    The equation is: F(x, y, z) = 0.1*x + 0.2*y - 0.1*z + 1.0
    """
    return 0.1*x + 0.2*y - 0.1*z + 1.0

def G(y1, x, y2):
    """
    Sample function G: R*R*R -> R.
    The fixed point condition is y = G(y, x, y).
    The equation is: G(y1, x, y2) = 0.1*y1 + 0.1*y2 + 0.3*x + 2.0.
    So, G(y, x, y) = 0.2*y + 0.3*x + 2.0
    """
    return 0.1*y1 + 0.1*y2 + 0.3*x + 2.0

def H(z, y, x):
    """
    Sample function H: R*R*R -> R.
    This function is part of a system that satisfies contraction conditions.
    The equation is: H(z, y, x) = 0.1*z - 0.2*y + 0.3*x + 3.0
    """
    return 0.1*z - 0.2*y + 0.3*x + 3.0

def find_tripled_fixed_point():
    """
    Finds the FGH-tripled fixed point using the method of successive approximations.
    """
    # Initial guess for the fixed point (x, y, z)
    x, y, z = 0.0, 0.0, 0.0

    # Number of iterations for convergence
    # For a contraction mapping, this iterative process is guaranteed to converge.
    num_iterations = 100

    for i in range(num_iterations):
        # Store the current values to use in this iteration's calculations
        x_prev, y_prev, z_prev = x, y, z
        
        # Calculate the next approximation for (x, y, z)
        x = F(x_prev, y_prev, z_prev)
        y = G(y_prev, x_prev, y_prev)
        z = H(z_prev, y_prev, x_prev)
    
    print("Found the tripled fixed point (x, y, z) after {} iterations:".format(num_iterations))
    print(f"({x:.4f}, {y:.4f}, {z:.4f})")
    print("-" * 50)
    
    print("Verifying the fixed point by substituting the values back into the equations:")
    
    # Equation 1: x = F(x, y, z)
    print("\n1. Condition x = F(x, y, z):")
    rhs_F = F(x, y, z)
    print(f"   {x:.4f} = 0.1 * {x:.4f} + 0.2 * {y:.4f} - 0.1 * {z:.4f} + 1.0")
    print(f"   {x:.4f} = {rhs_F:.4f}  (Difference: {abs(x - rhs_F):.2e})")

    # Equation 2: y = G(y, x, y)
    print("\n2. Condition y = G(y, x, y):")
    rhs_G = G(y, x, y)
    print(f"   {y:.4f} = 0.2 * {y:.4f} + 0.3 * {x:.4f} + 2.0")
    print(f"   {y:.4f} = {rhs_G:.4f}  (Difference: {abs(y - rhs_G):.2e})")

    # Equation 3: z = H(z, y, x)
    print("\n3. Condition z = H(z, y, x):")
    rhs_H = H(z, y, x)
    print(f"   {z:.4f} = 0.1 * {z:.4f} - 0.2 * {y:.4f} + 0.3 * {x:.4f} + 3.0")
    print(f"   {z:.4f} = {rhs_H:.4f}  (Difference: {abs(z - rhs_H):.2e})")

if __name__ == '__main__':
    find_tripled_fixed_point()
