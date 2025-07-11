import math

def calculate_furthest_distance_on_cone():
    """
    This function explains and calculates the expression for the furthest distance
    from a point P on a specific cone's base.
    
    The logic is based on unrolling the cone into a semicircle of radius 'd'
    and finding the longest straight-line distance from a point on the arc to a
    point in the semicircle.
    """
    print("The problem is solved by unrolling the cone into a semicircle of radius 'd'.")
    print("We then find the longest possible straight-line path on this 2D surface.")
    print("This path runs from the top of the semicircle's arc to one of its corners.")
    
    # The calculation is based on the distance formula between (0, d) and (d, 0).
    # Distance = sqrt( (d-0)^2 + (0-d)^2 ) = d * sqrt(1^2 + (-1)^2)
    
    a = 1
    b = -1
    
    # We present the final equation with its numeric components.
    print("\nThe final distance D is calculated as follows:")
    print(f"D = d * sqrt( ({a})^2 + ({b})^2 )")
    
    a_sq = a**2
    b_sq = b**2
    sum_sq = a_sq + b_sq
    
    print(f"D = d * sqrt( {a_sq} + {b_sq} )")
    print(f"D = d * sqrt( {sum_sq} )")
    
    coefficient = math.sqrt(sum_sq)
    print(f"\nThus, the furthest distance is d * sqrt(2).")
    print(f"The numerical value of the coefficient is approximately {coefficient:.4f}.")

calculate_furthest_distance_on_cone()
