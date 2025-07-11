def check_fgh_tripled_fixed_point_condition(k_f1, k_f2, k_f3, k_g1, k_g2, k_g3, k_h1, k_h2, k_h3):
    """
    Checks the sufficient condition for the existence and uniqueness of an FGH-tripled fixed point.

    The fixed point (x, y, z) is defined by the system:
    F(x, y, z) = x
    G(y, x, y) = y
    H(z, y, x) = z

    The condition is based on the Lipschitz constants for F, G, and H.

    Args:
        k_f1, k_f2, k_f3: Lipschitz constants for F w.r.t. its 1st, 2nd, and 3rd arguments.
        k_g1, k_g2, k_g3: Lipschitz constants for G w.r.t. its 1st, 2nd, and 3rd arguments.
        k_h1, k_h2, k_h3: Lipschitz constants for H w.r.t. its 1st, 2nd, and 3rd arguments.
    """
    print("--- FGH-Tripled Fixed Point Condition Check ---")
    print("The condition for a unique fixed point is: max(S_x, S_y, S_z) < 1")
    print("Where S_x, S_y, S_z are the sums of the Lipschitz constants affecting each space.")
    print("-" * 50)

    # Calculate the sums for each component of the product space
    sum_x = k_f1 + k_g2 + k_h3
    sum_y = k_f2 + k_g1 + k_g3 + k_h2
    sum_z = k_f3 + k_h1

    # Print the equation with the numbers substituted in
    print(f"S_x = k_f1 + k_g2 + k_h3 = {k_f1} + {k_g2} + {k_h3} = {sum_x:.4f}")
    print(f"S_y = k_f2 + k_g1 + k_g3 + k_h2 = {k_f2} + {k_g1} + {k_g3} + {k_h2} = {sum_y:.4f}")
    print(f"S_z = k_f3 + k_h1 = {k_f3} + {k_h1} = {sum_z:.4f}")
    print("-" * 50)

    # Final condition check
    overall_k = max(sum_x, sum_y, sum_z)
    
    print("Final Equation Check:")
    print(f"max({sum_x:.4f}, {sum_y:.4f}, {sum_z:.4f}) < 1")
    
    is_met = overall_k < 1

    print(f"Result: {overall_k:.4f} < 1 is {is_met}")

    if is_met:
        print("\nConclusion: The condition is MET. A unique FGH-tripled fixed point exists.")
    else:
        print("\nConclusion: The condition is NOT MET. This theorem does not guarantee a fixed point.")

# --- Example Usage ---
# Case 1: The condition is met
print("Example 1: A set of constants that satisfy the condition.")
check_fgh_tripled_fixed_point_condition(
    k_f1=0.2, k_f2=0.1, k_f3=0.1,
    k_g1=0.1, k_g2=0.2, k_g3=0.1,
    k_h1=0.1, k_h2=0.1, k_h3=0.2
)

print("\n" + "="*50 + "\n")

# Case 2: The condition is not met
print("Example 2: A set of constants that DO NOT satisfy the condition.")
check_fgh_tripled_fixed_point_condition(
    k_f1=0.4, k_f2=0.3, k_f3=0.2,
    k_g1=0.2, k_g2=0.5, k_g3=0.2,
    k_h1=0.3, k_h2=0.1, k_h3=0.4
)
