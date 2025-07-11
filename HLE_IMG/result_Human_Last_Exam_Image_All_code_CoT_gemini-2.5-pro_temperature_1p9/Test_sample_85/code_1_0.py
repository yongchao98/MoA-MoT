import math

def solve_cone_distance():
    """
    Calculates the furthest distance one can travel from a point P on the
    circumference of the base of the given cone.

    The method involves unwrapping the cone's surface into a semicircle and
    finding the maximum distance within that shape.
    """
    
    # The final expression for the distance is d * sqrt(2).
    # We need to output the numerical part of this expression as requested.
    coefficient = math.sqrt(2)

    # Print the explanation of the result.
    print("Based on the cone's properties, its unwrapped surface is a semicircle of radius 'd'.")
    print("The starting point P is on the arc of the semicircle.")
    print("The point on the cone's surface furthest from P corresponds to the point in the semicircle geometry that is furthest from P.")
    print("This furthest point lies on the diameter of the semicircle.")
    
    # Print the calculation steps.
    print("\nThe calculation for the squared distance is: d^2 + d^2 = 2 * d^2")
    print("Taking the square root gives the final distance.")
    
    # Print the final result and the number in the equation.
    print("\nThe final equation for the furthest distance is: d * sqrt(2)")
    print("The numerical value in this equation is sqrt(2).")
    print(f"sqrt(2) = {coefficient}")

solve_cone_distance()
