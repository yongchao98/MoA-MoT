import math

def solve_triangle_distance():
    """
    Calculates the shortest distance between the point with the maximum number
    and the point with the minimum number on the triangle grid.
    
    The user can modify the values of a, b, and c below.
    """
    
    # --- User-configurable values ---
    # The numbers placed on points A, B, C respectively.
    a = 10
    b = 5
    c = 10
    
    # The number of divisions on each side.
    # Note: The value of n does not affect the final distance.
    n = 100
    # ------------------------------------

    # The shortest distance depends on how many of the vertices A, B, C
    # share the maximum or minimum values. We can determine this by counting
    # the number of distinct values among a, b, and c.
    
    distinct_values = set([a, b, c])
    num_distinct = len(distinct_values)

    distance = 0.0

    if num_distinct == 1:
        # Case 1: a = b = c.
        # The maximum and minimum values in the grid are the same.
        # The shortest distance between a point with the max value and a
        # point with the min value is 0, as they can be the same point.
        distance = 0.0
        
    elif num_distinct == 2:
        # Case 2: Exactly two of a, b, c are equal.
        # This means either one vertex has an extremal value (e.g., a > b=c),
        # or two vertices share an extremal value (e.g., a=b > c).
        # In either scenario, one set of extremal points is a single vertex,
        # and the other set lies on the opposite edge. The shortest distance
        # is the altitude of the equilateral triangle of side length 1.
        # Altitude h = sqrt(side^2 - (side/2)^2) = sqrt(1 - 1/4) = sqrt(3/4).
        distance = math.sqrt(3) / 2.0
        
    elif num_distinct == 3:
        # Case 3: a, b, and c are all different.
        # The maximum value exists at only one vertex, and the minimum value
        # exists at only one other vertex. The distance between any two
        # distinct vertices of the equilateral triangle is the side length.
        distance = 1.0

    print(f"The values at vertices A, B, C are a={a}, b={b}, c={c}.")
    print(f"The shortest distance between the max and min value points is: {distance}")

solve_triangle_distance()