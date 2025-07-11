import math

def solve_cone_distance():
    """
    Calculates the furthest distance from a point P on the base of a cone,
    traveling only on the surface.
    """
    
    print("Step 1: Determine the properties of the cone from the problem statement.")
    print("Base diameter = d")
    print("Base radius, r = d / 2")
    print("Slant height, l = d (distance from apex to P)")
    
    print("\nStep 2: Unroll the cone's lateral surface into a sector of a circle.")
    print("The radius of the sector is the slant height, l = d.")
    print("The arc length of the sector is the base circumference, C = 2 * pi * r = pi * d.")
    
    print("\nStep 3: Calculate the angle of the sector in radians.")
    print("Sector angle, theta = Arc Length / Radius = (pi * d) / d = pi radians.")
    print("This means the unrolled surface is a semicircle of radius d.")
    
    print("\nStep 4: Find the furthest point on the unrolled surface.")
    print("Point P is on the arc of the semicircle. The furthest point is diametrically opposite on the arc.")
    print("The distance is the straight line connecting these two points, which is the diameter of the full circle.")
    
    print("\nStep 5: Calculate the final distance using the Law of Cosines for the triangle formed by the apex (V) and the two points (P, Q).")
    print("In the unrolled semicircle, V is the center, P and Q are endpoints of the diameter on the arc.")
    print("VP = d, VQ = d, and the angle PVQ = pi radians.")
    
    # Symbolic representation for the final equation printout
    d_str = 'd'
    val_2 = 2
    val_minus_1 = -1
    val_4 = 4
    final_dist = '2d'
    
    print("\nFinal Answer Equation:")
    print(f"Distance^2 = {d_str}^2 + {d_str}^2 - {val_2}*{d_str}*{d_str}*cos(pi)")
    print(f"Distance^2 = {d_str}^2 + {d_str}^2 - {val_2}*{d_str}^2*({val_minus_1})")
    print(f"Distance^2 = {val_2}{d_str}^2 + {val_2}{d_str}^2")
    print(f"Distance^2 = {val_4}{d_str}^2")
    print(f"Distance = sqrt({val_4}{d_str}^2) = {final_dist}")

solve_cone_distance()
<<<2d>>>