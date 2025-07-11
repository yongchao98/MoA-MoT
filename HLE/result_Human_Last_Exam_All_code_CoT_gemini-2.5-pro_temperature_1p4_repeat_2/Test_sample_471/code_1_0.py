def solve_critical_points_on_torus():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse theory.
    """

    # The Betti numbers for the 2-torus (T^2)
    # b_0: number of connected components
    b_0 = 1
    # b_1: number of 1D holes (non-contractible loops)
    b_1 = 2
    # b_2: number of 2D voids/cavities
    b_2 = 1
    
    betti_numbers = [b_0, b_1, b_2]

    # According to a result from Morse theory, the number of critical points of a smooth
    # function on a compact manifold is at least the sum of its Betti numbers.
    min_critical_points = sum(betti_numbers)

    # We demonstrate that this minimum is achievable.
    # Consider the height function on a torus standing vertically. It has:
    # c_0 = 1 minimum (at the bottom)
    # c_1 = 2 saddle points (on the inner rim)
    # c_2 = 1 maximum (at the top)
    # The total number of critical points for this function is 1 + 2 + 1 = 4.
    
    print("This problem can be solved using Morse Theory.")
    print("The minimal number of critical points for a smooth function on a manifold M")
    print("is greater than or equal to the sum of its Betti numbers (b_k).")
    print("\nFor the 2-torus (T^2):")
    print(f"Betti number b_0 (connected components) = {b_0}")
    print(f"Betti number b_1 (1D holes) = {b_1}")
    print(f"Betti number b_2 (2D voids) = {b_2}")
    
    print("\nThe minimal number of critical points is the sum:")
    # Using f-string to format the equation nicely.
    print(f"Minimal number = b_0 + b_1 + b_2 = {b_0} + {b_1} + {b_2} = {min_critical_points}")
    print("\nSince a function exists with exactly this many critical points (e.g., the height function on a standard torus), this is the minimum.")


solve_critical_points_on_torus()
<<<4>>>