def check_fgh_tripled_fixed_point_conditions(k_Fx, k_Fy, k_Fz, k_Gx, k_Gy, k_Gy_alt, k_Hx, k_Hy, k_Hz):
    """
    Checks the conditions for the existence and uniqueness of an FGH-tripled fixed point.

    The FGH-tripled fixed point is a solution (x, y, z) to the system:
        F(x, y, z) = x
        G(y, x, y) = y
        H(z, y, x) = z
    
    This function checks if the Lipschitz constants of F, G, and H satisfy the
    sufficient conditions derived from the Banach Contraction Mapping Principle.

    Args:
        k_Fx, k_Fy, k_Fz: Lipschitz constants for F w.r.t. x, y, z.
        k_Gx, k_Gy, k_Gy_alt: Lipschitz constants for G w.r.t. its 1st, 2nd, and 3rd arguments.
        k_Hx, k_Hy, k_Hz: Lipschitz constants for H w.r.t. x, y, z arguments.
    """
    print("--- Checking Conditions for a Unique FGH-Tripled Fixed Point ---\n")

    # Condition 1 (sum of coeffs for x-dimension)
    sum_x = k_Fx + k_Gx + k_Hx
    cond1_met = sum_x < 1
    print("Condition 1 (for x-dimension): k_Fx + k_Gx + k_Hx < 1")
    print(f"  Calculation: {k_Fx} + {k_Gx} + {k_Hx} = {sum_x}")
    print(f"  Result: {sum_x} < 1 is {cond1_met}\n")

    # Condition 2 (sum of coeffs for y-dimension)
    sum_y = k_Fy + k_Gy + k_Gy_alt + k_Hy
    cond2_met = sum_y < 1
    print("Condition 2 (for y-dimension): k_Fy + k_Gy + k_Gy_alt + k_Hy < 1")
    print(f"  Calculation: {k_Fy} + {k_Gy} + {k_Gy_alt} + {k_Hy} = {sum_y}")
    print(f"  Result: {sum_y} < 1 is {cond2_met}\n")

    # Condition 3 (sum of coeffs for z-dimension)
    sum_z = k_Fz + k_Hz
    cond3_met = sum_z < 1
    print("Condition 3 (for z-dimension): k_Fz + k_Hz < 1")
    print(f"  Calculation: {k_Fz} + {k_Hz} = {sum_z}")
    print(f"  Result: {sum_z} < 1 is {cond3_met}\n")

    # Final Conclusion
    print("--- Conclusion ---")
    if cond1_met and cond2_met and cond3_met:
        print("All conditions are met. A unique FGH-tripled fixed point is guaranteed to exist.")
    else:
        print("Not all conditions are met. A unique FGH-tripled fixed point is not guaranteed by this theorem.")


if __name__ == '__main__':
    # --- Example Usage ---
    # Define a set of Lipschitz constants for F, G, H.
    # Let's test a case where the conditions are met.
    print("Example 1: Constants that satisfy the conditions")
    check_fgh_tripled_fixed_point_conditions(
        k_Fx=0.3, k_Fy=0.1, k_Fz=0.4,
        k_Gx=0.2, k_Gy=0.2, k_Gy_alt=0.1,
        k_Hx=0.2, k_Hy=0.2, k_Hz=0.4
    )
    print("\n" + "="*60 + "\n")
    
    # Let's test a case where the conditions are NOT met.
    print("Example 2: Constants that do not satisfy the conditions")
    check_fgh_tripled_fixed_point_conditions(
        k_Fx=0.5, k_Fy=0.3, k_Fz=0.6,
        k_Gx=0.4, k_Gy=0.2, k_Gy_alt=0.3,
        k_Hx=0.2, k_Hy=0.3, k_Hz=0.5
    )