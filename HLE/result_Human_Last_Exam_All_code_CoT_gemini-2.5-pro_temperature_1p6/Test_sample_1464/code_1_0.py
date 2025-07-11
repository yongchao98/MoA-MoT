import math

def solve_and_present():
    """
    This function solves the problem by following the analytical steps
    of factoring the quartic polynomial, solving the resulting quadratics,
    and sorting the final roots.
    """
    print("The polynomial is factored into two quadratic equations.\n")

    # First quadratic equation: X^2 - s1*X + p1 = 0
    s1_val = math.sqrt(14) + 2 * math.sqrt(11)
    p1_val = 2 * math.sqrt(154)
    s1_str = "sqrt(14) + 2*sqrt(11)"
    p1_str = "2*sqrt(154)"
    print(f"Equation 1: X^2 - ({s1_str})*X + ({p1_str}) = 0")
    
    # Discriminant delta1 = s1^2 - 4*p1
    # s1^2 = 14 + 44 + 4*sqrt(154) = 58 + 4*sqrt(154)
    # 4*p1 = 8*sqrt(154)
    # delta1 = 58 - 4*sqrt(154)
    # sqrt(delta1) simplifies to 2*sqrt(11) - sqrt(14)
    sqrt_delta1_val = 2 * math.sqrt(11) - math.sqrt(14)
    
    # Roots for the first quadratic
    root1_val = (s1_val + sqrt_delta1_val) / 2
    root2_val = (s1_val - sqrt_delta1_val) / 2
    
    root1_str = "2*sqrt(11)"
    root2_str = "sqrt(14)"
    print(f"The roots are {root1_str} and {root2_str}.\n")
    
    # Second quadratic equation: X^2 - s2*X + p2 = 0
    s2_val = math.sqrt(34) + 2 * math.sqrt(6)
    p2_val = 4 * math.sqrt(51)
    s2_str = "sqrt(34) + 2*sqrt(6)"
    p2_str = "4*sqrt(51)"
    print(f"Equation 2: X^2 - ({s2_str})*X + ({p2_str}) = 0")
    
    # Discriminant delta2 = s2^2 - 4*p2
    # s2^2 = 34 + 24 + 4*sqrt(204) = 58 + 8*sqrt(51)
    # 4*p2 = 16*sqrt(51)
    # delta2 = 58 - 8*sqrt(51)
    # sqrt(delta2) simplifies to sqrt(34) - 2*sqrt(6)
    sqrt_delta2_val = math.sqrt(34) - 2 * math.sqrt(6)

    # Roots for the second quadratic
    root3_val = (s2_val + sqrt_delta2_val) / 2
    root4_val = (s2_val - sqrt_delta2_val) / 2

    root3_str = "sqrt(34)"
    root4_str = "2*sqrt(6)"
    print(f"The roots are {root3_str} and {root4_str}.\n")

    # Now, we collect all symbolic roots and their numerical values to sort them
    all_roots = {
        root1_str: root1_val,
        root2_str: root2_val,
        root3_str: root3_val,
        root4_str: root4_val,
    }

    # Sort the roots based on their numerical value (the dict's value)
    sorted_roots = sorted(all_roots.items(), key=lambda item: item[1])

    print("The four roots of the original equation in increasing order are:")
    # Print the sorted symbolic representations
    final_root_strings = [item[0] for item in sorted_roots]
    print(", ".join(final_root_strings))

solve_and_present()