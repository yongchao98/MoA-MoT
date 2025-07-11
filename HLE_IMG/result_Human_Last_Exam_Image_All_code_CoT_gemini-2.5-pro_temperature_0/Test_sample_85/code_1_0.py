import math

def solve_cone_distance():
    """
    This function explains the steps to find the furthest distance on a cone's surface
    and prints the final expression for the answer.
    """
    print("To find the furthest distance from point P on the cone's surface, we unroll the cone.")
    
    # Step 1: Describe the unrolled shape
    print("\nStep 1: Determine the shape of the unrolled cone surface.")
    print("The cone has a base diameter 'd' and a slant height 'l' also equal to 'd'.")
    print("When unrolled, the cone's lateral surface forms a sector of a circle.")
    print("The sector's radius is the slant height 'd'.")
    print("The sector's arc length is the base circumference 'pi * d'.")
    print("The sector's angle is (pi * d) / d = pi radians (180 degrees).")
    print("Result: The unrolled surface is a semicircle of radius 'd'.")

    # Step 2: Identify the furthest point
    print("\nStep 2: Locate the furthest point on the semicircle.")
    print("Let the apex be the center of the semicircle. Point P corresponds to the two ends of the diameter.")
    print("The point 'Q' on the surface furthest from P must be equidistant from these two ends.")
    print("This places Q at the peak of the semicircle's arc.")

    # Step 3: Calculate the distance
    print("\nStep 3: Calculate the distance.")
    print("The distance is the hypotenuse of a right-angled triangle where the two legs are the radius 'd'.")
    print("Using the Pythagorean theorem: Distance = sqrt(d^2 + d^2) = sqrt(2 * d^2).")

    # Final Answer
    print("\n--------------------------------------------------")
    print("The final equation for the furthest distance is:")
    
    # The variable is 'd'. The number in the equation is 2.
    variable_name = 'd'
    number_in_equation = 2
    
    print(f"{variable_name} * sqrt({number_in_equation})")
    print("--------------------------------------------------")

solve_cone_distance()