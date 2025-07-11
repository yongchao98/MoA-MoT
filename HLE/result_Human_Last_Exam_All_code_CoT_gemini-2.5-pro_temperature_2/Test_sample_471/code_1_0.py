def solve_critical_points_on_torus():
    """
    Calculates the minimal number of critical points for a smooth function on a 2-torus
    using Morse theory and Betti numbers.
    """
    # The Betti numbers for the 2-torus (T^2) describe its topological properties.
    # b_0 is the number of connected components.
    b_0 = 1
    # b_1 is the number of 1-dimensional "holes".
    b_1 = 2
    # b_2 is the number of 2-dimensional "voids".
    b_2 = 1

    print("To find the minimal number of critical points for a smooth function f: T^2 -> R, we use Morse Theory.")
    print("The total number of critical points is bounded below by the sum of the Betti numbers of the manifold T^2.\n")
    print(f"The Betti numbers for the 2-torus are:")
    print(f"b_0 = {b_0} (number of connected components)")
    print(f"b_1 = {b_1} (number of 1-dimensional holes)")
    print(f"b_2 = {b_2} (number of 2-dimensional voids)\n")
    
    # According to the Morse inequalities, the total number of critical points must be at least the sum of the Betti numbers.
    min_critical_points = b_0 + b_1 + b_2
    
    print("The minimal number of critical points is the sum of these Betti numbers.")
    print("Final Calculation:")
    print(f"{b_0} + {b_1} + {b_2} = {min_critical_points}")
    print("\nThis bound is achievable. For instance, a standard height function on a torus can be constructed to have exactly 1 minimum, 2 saddle points, and 1 maximum, totaling 4 critical points.")
    print(f"Thus, the minimal number of critical points is {min_critical_points}.")

solve_critical_points_on_torus()