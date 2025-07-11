def solve_critical_points_on_torus():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse theory.
    """
    # The Betti numbers for the 2-torus (T^2)
    # b_0: number of connected components
    b0 = 1
    # b_1: number of "holes" (one-dimensional cycles)
    b1 = 2
    # b_2: number of "voids" (two-dimensional cycles)
    b2 = 1

    # According to the Morse inequalities, the number of critical points c_k of
    # index k must satisfy c_k >= b_k for a Morse function.
    # The total number of critical points N must therefore be at least the
    # sum of the Betti numbers. This lower bound is the actual minimum.
    min_critical_points = b0 + b1 + b2

    print("The Betti numbers of the 2-torus are:")
    print(f"b_0 = {b0}")
    print(f"b_1 = {b1}")
    print(f"b_2 = {b2}")
    print("\nThe minimal number of critical points is the sum of these Betti numbers:")
    print(f"{b0} + {b1} + {b2} = {min_critical_points}")
    print("\nThis corresponds to 1 minimum (index 0), 2 saddle points (index 1), and 1 maximum (index 2).")


solve_critical_points_on_torus()