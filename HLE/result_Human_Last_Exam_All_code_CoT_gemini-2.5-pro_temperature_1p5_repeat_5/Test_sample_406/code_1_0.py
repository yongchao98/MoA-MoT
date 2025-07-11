import numpy as np

def find_tripled_fixed_point():
    """
    This function demonstrates the conditions for an FGH-tripled fixed point.
    
    It defines a specific system of functions F, G, H over the real numbers 
    (which is a complete metric space) that satisfy the contractivity conditions.
    Then, it calculates the unique tripled fixed point for this system.
    """
    
    # 1. Define the contractive mappings and coefficients.
    # Let X = Y = Z = R, with the metric d(a,b) = |a-b|.
    # We define linear functions which automatically satisfy the contractive inequalities.
    # F(x, y, z) = a1*x + a2*y + a3*z + const_f
    # G(y, x, y) = b1*y + b2*x + const_g
    # H(z, y, x) = c1*z + c2*y + c3*x + const_h
    
    # Choose coefficients
    a1, a2, a3 = 0.1, 0.1, 0.1
    b1, b2 = 0.2, 0.1
    c1, c2, c3 = 0.1, 0.2, 0.1
    
    # Choose constant terms
    const_f, const_g, const_h = 1.0, 2.0, 3.0
    
    # 2. Check the contraction factor condition.
    k_x = a1 + b2 + c3
    k_y = a2 + b1 + c2
    k_z = a3 + c1
    k = max(k_x, k_y, k_z)
    
    print("--- Conditions Check ---")
    print(f"Coefficients for F: a1={a1}, a2={a2}, a3={a3}")
    print(f"Coefficients for G: b1={b1}, b2={b2}")
    print(f"Coefficients for H: c1={c1}, c2={c2}, c3={c3}")
    print(f"Contraction factor k = max({k_x:.2f}, {k_y:.2f}, {k_z:.2f}) = {k:.2f}")

    if k >= 1:
        print("Condition not met: k >= 1. A unique fixed point is not guaranteed.")
        return
    else:
        print("Condition met: k < 1. A unique fixed point exists.")

    # 3. Find the fixed point (x, y, z).
    # The fixed point equations are:
    # x = a1*x + a2*y + a3*z + const_f  => (1-a1)x - a2*y - a3*z = const_f
    # y = b1*y + b2*x + const_g          => -b2*x + (1-b1)y       = const_g
    # z = c1*z + c2*y + c3*x + const_h  => -c3*x - c2*y + (1-c1)z = const_h
    
    # This is a linear system Ax = B
    A = np.array([
        [1 - a1, -a2, -a3],
        [-b2, 1 - b1, 0],
        [-c3, -c2, 1 - c1]
    ])
    
    B = np.array([const_f, const_g, const_h])
    
    # Solve for the fixed point (x, y, z)
    try:
        x_fp, y_fp, z_fp = np.linalg.solve(A, B)
        print(f"\nThe unique tripled fixed point (x, y, z) is: ({x_fp:.4f}, {y_fp:.4f}, {z_fp:.4f})")
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # 4. Verify the solution by plugging the numbers into the equations.
    print("\n--- Verification of the Final Equations ---")
    
    # Equation for F(x, y, z) = x
    lhs_f = a1*x_fp + a2*y_fp + a3*z_fp + const_f
    print(f"F({x_fp:.3f}, {y_fp:.3f}, {z_fp:.3f}) = {a1}*{x_fp:.3f} + {a2}*{y_fp:.3f} + {a3}*{z_fp:.3f} + {const_f} = {lhs_f:.3f}")
    print(f"Result: {lhs_f:.3f}  -->  x = {x_fp:.3f}\n")
    
    # Equation for G(y, x, y) = y
    lhs_g = b1*y_fp + b2*x_fp + const_g
    print(f"G({y_fp:.3f}, {x_fp:.3f}, {y_fp:.3f}) = {b1}*{y_fp:.3f} + {b2}*{x_fp:.3f} + {const_g} = {lhs_g:.3f}")
    print(f"Result: {lhs_g:.3f}  -->  y = {y_fp:.3f}\n")
    
    # Equation for H(z, y, x) = z
    lhs_h = c1*z_fp + c2*y_fp + c3*x_fp + const_h
    print(f"H({z_fp:.3f}, {y_fp:.3f}, {x_fp:.3f}) = {c1}*{z_fp:.3f} + {c2}*{y_fp:.3f} + {c3}*{x_fp:.3f} + {const_h} = {lhs_h:.3f}")
    print(f"Result: {lhs_h:.3f}  -->  z = {z_fp:.3f}")


if __name__ == '__main__':
    find_tripled_fixed_point()
