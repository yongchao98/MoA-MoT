import math

def solve_tangent_circles():
    """
    Solves for r^2 based on the geometric constraints of two tangent circles.
    """
    # Step 1: Analyze the geometry of the lines
    print("Step 1: Analyze the two lines y = x + 1 and y = -x + 5.")
    print("The lines can be written in the standard form Ax + By + C = 0 as:")
    print("L1: x - y + 1 = 0")
    print("L2: x + y - 5 = 0")
    print("The slope of L1 is 1 and the slope of L2 is -1. Since 1 * -1 = -1, the lines are perpendicular.")
    print("Their intersection is found by solving the system of equations, which gives the point (2, 3).")
    print("The angle bisectors of two perpendicular lines are x=2 and y=3.")
    print("The centers of the circles C and D must lie on one of these bisectors.")
    print("-" * 30)

    # Step 2 & 3: Relate center to radius on a bisector
    print("Step 2: Relate a circle's center to its radius.")
    print("Let's assume the centers lie on the bisector y = 3. Let a center be at (x, 3).")
    print("The radius is the distance from the center to either line.")
    print("Using the distance formula from a point (x0, y0) to a line Ax+By+C=0: d = |A*x0 + B*y0 + C| / sqrt(A^2 + B^2).")
    print("For L2: radius = |1*x + 1*3 - 5| / sqrt(1^2 + 1^2) = |x - 2| / sqrt(2).")
    print("-" * 30)

    # Step 4: Use information about the second circle (D)
    print("Step 3: Find the position of the center D of the second circle.")
    print("The second circle has radius R = 2. Its center is D = (xd, 3).")
    print("Using the radius formula: 2 = |xd - 2| / sqrt(2).")
    print("This gives |xd - 2| = 2 * sqrt(2).")
    print("Let's place D at xd = 2 + 2 * sqrt(2). The choice of the other solution, xd = 2 - 2*sqrt(2), would lead to a symmetric result.")
    print("-" * 30)

    # Step 5 & 6: Apply tangency and solve for r
    print("Step 4: Use the tangency condition between the two circles to find r.")
    print("The first circle has center C = (xc, 3) and radius r = |xc - 2| / sqrt(2).")
    print("The circles are tangent, so the distance between centers equals the sum of radii: distance(C, D) = r + R.")
    print("|xc - xd| = r + 2")
    print("|xc - (2 + 2*sqrt(2))| = |xc - 2| / sqrt(2) + 2")
    print("\nThe problem refers to a 'first circle' (C) and a 'second circle' (D).")
    print("This implies an ordering. We assume the first circle C is closer to the line intersection point (2, 3) than D is.")
    print("This means |xc - 2| < |xd - 2|, so |xc - 2| < 2*sqrt(2).")
    print("This implies xc is between 2 - 2*sqrt(2) and 2 + 2*sqrt(2).")
    print("With this condition, the equation becomes:")
    print("(2 + 2*sqrt(2)) - xc = (xc - 2) / sqrt(2) + 2")
    print("Solving for xc: 2*sqrt(2) - xc = (xc-2)/sqrt(2)")
    print("2*2 - xc*sqrt(2) = xc - 2")
    print("4 + 2 = xc * (1 + sqrt(2)) => 6 = xc * (1 + sqrt(2))")
    print("xc = 6 / (1 + sqrt(2)) = 6 * (sqrt(2) - 1) = 6*sqrt(2) - 6")
    print("\nNow we find the radius r using xc:")
    print("r = (xc - 2) / sqrt(2) because C is to the right of x=2 (since 6*sqrt(2)-6 > 2).")
    print("r = ( (6*sqrt(2) - 6) - 2 ) / sqrt(2)")
    print("r = (6*sqrt(2) - 8) / sqrt(2)")
    print("r = 6 - 8/sqrt(2) = 6 - 4*sqrt(2)")
    print("-" * 30)

    # Step 7: Calculate r^2
    print("Step 5: Calculate the final value of r^2.")
    print("r^2 = (6 - 4*sqrt(2))^2")
    r_squared_term1 = 6**2
    r_squared_term2 = 2 * 6 * 4
    r_squared_term3 = (4**2) * 2
    print(f"r^2 = {r_squared_term1} - {r_squared_term2}*sqrt(2) + {r_squared_term3}")
    final_a = r_squared_term1 + r_squared_term3
    final_b = r_squared_term2
    print(f"r^2 = {final_a} - {final_b}*sqrt(2)")
    
solve_tangent_circles()