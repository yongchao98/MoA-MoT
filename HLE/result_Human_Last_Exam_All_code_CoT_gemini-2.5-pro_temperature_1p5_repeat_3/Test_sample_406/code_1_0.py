def check_fgh_tripled_fixed_point_conditions(coeffs):
    """
    Checks if the given contraction coefficients guarantee a unique FGH-tripled fixed point.

    Args:
        coeffs (dict): A dictionary containing the contraction coefficients.
            Keys: 'f_x', 'f_y', 'f_z' for function F(x, y, z)
                  'g_x', 'g_y1', 'g_y2' for function G(y, x, y)
                  'h_x', 'h_y', 'h_z' for function H(z, y, x)
    """
    f_x, f_y, f_z = coeffs['f_x'], coeffs['f_y'], coeffs['f_z']
    g_x, g_y1, g_y2 = coeffs['g_x'], coeffs['g_y1'], coeffs['g_y2']
    h_x, h_y, h_z = coeffs['h_x'], coeffs['h_y'], coeffs['h_z']

    print("Analyzing the conditions for an FGH-tripled fixed point.")
    print("-" * 60)
    print("The fixed point (x, y, z) is defined by the system:")
    print("  F(x, y, z) = x")
    print("  G(y, x, y) = y")
    print("  H(z, y, x) = z")
    print("-" * 60)
    print("Based on the Banach Fixed-Point Theorem on the product space X*Y*Z,")
    print("a unique fixed point is guaranteed if the following inequalities hold:")

    # Calculate the contraction factors for each dimension
    k_x = f_x + g_x + h_x
    k_y = f_y + g_y1 + g_y2 + h_y
    k_z = f_z + h_z
    
    # Check and print the condition for X
    print("\n1. Condition on X component:")
    print(f"   Equation: f_x + g_x + h_x < 1")
    print(f"   With your numbers: {f_x} + {g_x} + {h_x} = {k_x:.4f}")
    is_k_x_valid = k_x < 1
    print(f"   Result: {k_x:.4f} < 1 is {is_k_x_valid}")

    # Check and print the condition for Y
    print("\n2. Condition on Y component:")
    print(f"   Equation: f_y + g_y1 + g_y2 + h_y < 1")
    print(f"   With your numbers: {f_y} + {g_y1} + {g_y2} + {h_y} = {k_y:.4f}")
    is_k_y_valid = k_y < 1
    print(f"   Result: {k_y:.4f} < 1 is {is_k_y_valid}")
    
    # Check and print the condition for Z
    print("\n3. Condition on Z component:")
    print(f"   Equation: f_z + h_z < 1")
    print(f"   With your numbers: {f_z} + {h_z} = {k_z:.4f}")
    is_k_z_valid = k_z < 1
    print(f"   Result: {k_z:.4f} < 1 is {is_k_z_valid}")
    
    # Final conclusion
    print("-" * 60)
    if is_k_x_valid and is_k_y_valid and is_k_z_valid:
        print("Conclusion: All conditions are met. A unique FGH-tripled fixed point is GUARANTEED.")
    else:
        print("Conclusion: At least one condition is not met. A unique FGH-tripled fixed point is NOT guaranteed by this theorem.")

if __name__ == '__main__':
    # Example 1: Conditions are met
    print("--- Example 1: A set of coefficients that works ---")
    coefficients_success = {
        'f_x': 0.2, 'f_y': 0.1, 'f_z': 0.3,
        'g_x': 0.3, 'g_y1': 0.2, 'g_y2': 0.1,
        'h_x': 0.2, 'h_y': 0.2, 'h_z': 0.4
    }
    check_fgh_tripled_fixed_point_conditions(coefficients_success)

    print("\n" + "="*60 + "\n")

    # Example 2: Conditions are not met
    print("--- Example 2: A set of coefficients that fails ---")
    coefficients_fail = {
        'f_x': 0.5, 'f_y': 0.2, 'f_z': 0.4,
        'g_x': 0.3, 'g_y1': 0.2, 'g_y2': 0.1,
        'h_x': 0.3, 'h_y': 0.3, 'h_z': 0.5
    }
    check_fgh_tripled_fixed_point_conditions(coefficients_fail)
