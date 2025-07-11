import math

def solve_circle_problem():
    """
    Solves for r^2 based on the geometric constraints of the two circles.
    """
    # Step 1: Analyze the lines
    # The lines are L1: y = x + 1 and L2: y = -x + 5.
    # Their slopes are 1 and -1, so they are perpendicular.
    # Their intersection point I is (2, 3).
    # The angle bisectors are the lines x = 2 and y = 3.

    # Step 2: Relate radius to center position
    # For a circle with radius 'r' and center 'C' tangent to both lines,
    # the distance from C to the intersection I(2,3) is d(C, I) = r * sqrt(2).

    # Step 3 & 4: Set up conditions for both circles
    # Let circle C have radius r and circle D have radius R.
    R = 2
    # We assume both circles are tangent to both lines.
    # d(C, I) = r * sqrt(2)
    # d(D, I) = R * sqrt(2) = 2 * sqrt(2)
    # For the circles to be tangent, their centers C and D must lie on the same bisector line.
    # The distance between their centers, d(C,D), depends on their relative positions.
    # Case 1: C and D are on the same side of I. d(C,D) = |d(C,I) - d(D,I)| = |r - R|*sqrt(2)
    # Case 2: C and D are on opposite sides of I. d(C,D) = d(C,I) + d(D,I) = (r + R)*sqrt(2)
    # The tangency condition between circles is d(C,D) = r + R (external) or d(C,D) = |r - R| (internal).

    # Step 5: Solve for r
    # We test the possible combinations. The only one that yields a unique integer solution is:
    # Case 1 with internal tangency:
    # |r - R| * sqrt(2) = |r - R|
    # This equation holds if |r - R| = 0, which means r = R.
    # The other possibilities either have no solution or result in multiple, non-integer solutions,
    # suggesting this is the intended interpretation.
    
    # Set r equal to R
    r = R

    # Calculate r-squared
    r_squared = r**2

    print("The problem implies both circles are tangent to the two lines.")
    print("This leads to a relationship between the distance between their centers (d) and their radii (r and R=2):")
    print("d = |r - 2| * sqrt(2)  or  d = (r + 2) * sqrt(2)")
    print("The tangency condition for the circles is d = |r - 2| or d = r + 2.")
    print("The only case that provides a unique solution is |r - 2| * sqrt(2) = |r - 2|, which implies r=2.")
    
    print("\nSolving for r, we find:")
    print(f"r = {r}")
    
    print("\nThe problem asks for the value of r^2. The final equation is:")
    print(f"{int(r)} * {int(r)} = {int(r_squared)}")

solve_circle_problem()