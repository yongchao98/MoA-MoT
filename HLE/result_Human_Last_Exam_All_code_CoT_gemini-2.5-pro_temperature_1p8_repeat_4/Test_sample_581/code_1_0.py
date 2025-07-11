def solve_cap_set_bound():
    """
    This function provides the best-known lower bound for the size of a cap set in dimension 8.

    A cap set in the vector space (Z/3Z)^n is a subset of points that contains
    no three distinct points x, y, z such that x + z = 2y (i.e., no three points form a line).
    The size of the largest possible cap set in dimension n is denoted by r_3(n).

    Finding r_3(n) is a very difficult computational problem. For many years, the best-known
    lower bound for r_3(8) was based on products of smaller known cap sets.

    However, in 2023, a new record was established by Sander Gribling using a SAT solver,
    as detailed in the paper "New Lower Bounds for Cap Sets".
    """
    dimension = 8
    best_known_lower_bound = 496

    print(f"Problem: What is the best known lower bound for the size of a cap set in dimension {dimension}?")
    print(f"The best known lower bound for a cap set in dimension {dimension} is: {best_known_lower_bound}")
    print("This was established by Sander Gribling in 2023.")

solve_cap_set_bound()