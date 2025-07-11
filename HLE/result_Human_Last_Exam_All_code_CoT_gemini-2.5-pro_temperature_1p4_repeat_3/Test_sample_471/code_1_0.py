def solve_minimal_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse inequalities.
    """
    # Betti numbers for the 2-torus (T^2)
    # b0: number of connected components
    # b1: number of 1-dimensional holes
    # b2: number of 2-dimensional holes
    b0 = 1
    b1 = 2
    b2 = 1

    print(f"The Betti numbers for the 2-torus are:")
    print(f"b0 = {b0}")
    print(f"b1 = {b1}")
    print(f"b2 = {b2}\n")

    # According to Morse theory, any compact manifold must have at least one
    # minimum and one maximum. For a Morse function, the number of critical points
    # of index k (c_k) is bounded by the Betti numbers (b_k).

    # Minimal number of minima (index 0)
    # c0 >= b0
    c0_min = b0

    # Minimal number of maxima (index 2)
    # This is equivalent to minima for the function -f.
    # c2 >= b2
    c2_min = b2

    # Minimal number of saddle points (index 1) from the second Morse inequality:
    # c1 - c0 >= b1 - b0  =>  c1 >= (b1 - b0) + c0
    # To find the absolute minimum c1, we use the minimum c0.
    c1_min = (b1 - b0) + c0_min
    
    print("From the Morse inequalities, we can derive the minimum number of critical points of each type:")
    print(f"Minimal minima (c0) >= b0 = {c0_min}")
    print(f"Minimal maxima (c2) >= b2 = {c2_min}")
    print(f"Minimal saddles (c1) >= (b1 - b0) + c0 = ({b1} - {b0}) + {c0_min} = {c1_min}\n")

    # The minimal total number of critical points is the sum of these minimums.
    total_min = c0_min + c1_min + c2_min

    print("The minimal total number of critical points is the sum of these minimums:")
    print(f"{c0_min} (minima) + {c1_min} (saddles) + {c2_min} (maxima) = {total_min}")

solve_minimal_critical_points()