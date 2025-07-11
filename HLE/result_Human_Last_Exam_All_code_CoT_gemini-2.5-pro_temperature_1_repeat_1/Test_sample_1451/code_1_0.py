def solve_and_print():
    """
    Calculates T(n) for the given values and prints the equations and final results.
    """
    # For n=2:
    n2 = 2
    # Since n=2 is even, the formula T(n) = 2n - 1 applies.
    t2 = 2 * n2 - 1
    print(f"For n = {n2} (even):")
    print(f"T({n2}) = 2 * {n2} - 1 = {t2}")
    print("-" * 20)

    # For n=3:
    n3 = 3
    # Since n=3 is odd, the formula T(n) = 2n - 2 applies.
    t3 = 2 * n3 - 2
    print(f"For n = {n3} (odd):")
    print(f"T({n3}) = 2 * {n3} - 2 = {t3}")
    print("-" * 20)

    # For n=1234:
    n1234 = 1234
    # Since n=1234 is even, the formula T(n) = 2n - 1 applies.
    t1234 = 2 * n1234 - 1
    print(f"For n = {n1234} (even):")
    print(f"T({n1234}) = 2 * {n1234} - 1 = {t1234}")
    print("-" * 20)

    # For n=6712:
    n6712 = 6712
    # Since n=6712 is even, the formula T(n) = 2n - 1 applies.
    t6712 = 2 * n6712 - 1
    print(f"For n = {n6712} (even):")
    print(f"T({n6712}) = 2 * {n6712} - 1 = {t6712}")
    print("-" * 20)

    # Print the final answer as a comma-separated list
    print("\nThe requested values for T(2), T(3), T(1234), and T(6712) are:")
    print(f"{t2},{t3},{t1234},{t6712}")

solve_and_print()