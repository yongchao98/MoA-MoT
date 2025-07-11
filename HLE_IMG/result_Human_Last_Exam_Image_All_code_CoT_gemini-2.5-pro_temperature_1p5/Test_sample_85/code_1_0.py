import math

def solve_cone_distance():
    """
    Calculates the furthest distance from a point P on the base of a specific cone.
    
    The problem is solved by unrolling the cone into a semicircle and finding the
    maximum distance from a point on the arc to any other point in the semicircle.
    """
    
    # The problem defines the parameters in terms of 'd'.
    # Base diameter = d
    # Slant height L = d
    # Base radius r = d / 2

    # Step 1: Calculate the angle of the unrolled sector.
    # Circumference of base C = 2 * pi * r = pi * d
    # Arc length of sector S = C = pi * d
    # Radius of sector R = L = d
    # Angle theta = S / R = (pi * d) / d = pi radians (180 degrees)
    # The unrolled shape is a semicircle of radius d.
    
    print("Step-by-step derivation:")
    print("1. The base radius is r = d/2 and the slant height is L = d.")
    print("2. Unrolling the cone gives a sector of a circle of radius L=d.")
    print("3. The arc length of this sector is the cone's base circumference, C = pi*d.")
    print("4. The angle of the sector is theta = (arc length)/radius = (pi*d)/d = pi radians.")
    print("5. This means the unrolled surface is a semicircle of radius d.")
    print("6. Let the starting point P be at the top of the semicircle's arc.")
    print("7. The furthest point Q is at the corner of the semicircle's diameter.")
    print("8. We calculate the distance using the Pythagorean theorem.")
    print("   Distance^2 = (radius)^2 + (radius)^2 = d^2 + d^2 = 2*d^2")
    print("   Distance = sqrt(2*d^2) = d * sqrt(2)")
    
    print("\nThe final answer is the expression for the distance.")
    
    # The final equation is Distance = d * sqrt(2).
    # The prompt asks to output each number in the final equation. The number here is 2.
    number_in_equation = 2
    
    print("\nFinal Equation:")
    print(f"d * sqrt({number_in_equation})")

solve_cone_distance()
