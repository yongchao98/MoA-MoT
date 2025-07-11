import sys

def solve():
    """
    Solves the geometry problem to find the maximum value of n.
    """
    
    # Let S = { P_1, P_2, ... P_n } be a set of n points equidistant from another point O.
    # This means the points of S lie on a circle with center O.
    # Let T = S U {O}.
    # We can draw 9 straight lines to connect these points.

    # Number of straight lines available.
    num_lines = 9

    # A straight line can intersect a circle at most this many times.
    max_points_per_line_on_circle = 2

    # The condition is that any point in T can be reached from any other point in T 
    # by travelling along at most 2 of the 9 lines.

    # To find the maximum value of n, we first establish an upper bound.
    # The n points of S are on a circle. Each of the 9 lines can contain at most 2 points of S.
    # Therefore, the total number of points 'n' cannot exceed the number of lines multiplied by 
    # the maximum number of points per line.
    
    # The equation for the upper bound of n:
    # n <= num_lines * max_points_per_line_on_circle
    max_n = num_lines * max_points_per_line_on_circle

    print("Step 1: Establishing the upper bound for n.")
    print("The 'n' points of set S all lie on a single circle.")
    print("A single straight line can intersect a circle at most twice.")
    print(f"With {num_lines} lines available, the maximum possible value for n is given by the equation:")
    print(f"n <= {num_lines} * {max_points_per_line_on_circle}")
    print(f"n <= {max_n}")
    print("-" * 20)

    print("Step 2: Proving the maximum value is achievable.")
    print("To show that n=18 is the maximum, we must demonstrate a valid configuration.")
    print("Consider a configuration where all 9 lines pass through the center point O.")
    print("Each of the 9 lines intersects the circle at 2 unique points.")
    print(f"This construction gives us n = {num_lines} * {max_points_per_line_on_circle} = {max_n} points in the set S.")
    print("\nChecking the connectivity condition for this configuration:")
    print("  - Path between O and any point P_i in S:")
    print("    P_i lies on one of the lines (say L_k). Since all lines pass through O, O is also on L_k. The path uses just 1 line.")
    print("  - Path between two points P_i and P_j in S:")
    print("    P_i lies on line L_a and P_j lies on line L_b. All lines intersect at point O.")
    print("    A path exists from P_i to O (on line L_a) and from O to P_j (on line L_b). The path uses at most 2 lines.")
    print("-" * 20)
    
    print("Conclusion:")
    print("The upper bound for n is 18, and a valid configuration exists for n = 18.")
    print(f"Therefore, the maximum value of n is {max_n}.")

solve()