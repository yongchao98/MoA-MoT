import math

def solve():
    """
    This function determines if the five legs of a chair can touch a spherical surface simultaneously
    by checking if the five coplanar points representing the leg tips are concyclic.
    """

    # Step 1: Define the points for the five leg tips based on the problem description.
    # The coordinates imply the leg tips are in a single plane.
    points = {
        'p1': (0, 0),
        'p2': (2, 0),
        'p3': (2, 2),
        'p4': (0, 2),
        'p5': (1, 4)
    }
    point_list = list(points.values())
    point_names = list(points.keys())

    print("The five leg tip positions are:")
    for name, p in points.items():
        print(f"{name}: {p}")
    print("-" * 40)

    # Step 2: Determine the equation of a circle from three non-collinear points.
    # The equation is (x - h)^2 + (y - k)^2 = r^2.
    # We can find the center (h, k) by finding the intersection of the perpendicular
    # bisectors of two chords. Let's use the chord between p1 and p2, and p1 and p4.
    p1 = points['p1']
    p2 = points['p2']
    p4 = points['p4']

    # Perpendicular bisector of p1(0,0) and p2(2,0) is the vertical line x = (0+2)/2 = 1.
    h = (p1[0] + p2[0]) / 2.0
    # Perpendicular bisector of p1(0,0) and p4(0,2) is the horizontal line y = (0+2)/2 = 1.
    k = (p1[1] + p4[1]) / 2.0
    center = (h, k)

    # Step 3: Calculate the radius squared (r^2) using one of the points (e.g., p1).
    r_squared = (p1[0] - h)**2 + (p1[1] - k)**2

    print("To check if the points are concyclic, we define a circle using p1, p2, and p4.")
    print(f"The equation of the circle is (x - {h})^2 + (y - {k})^2 = {r_squared}")
    print("-" * 40)

    # Step 4: Check if all five points lie on this circle.
    print("Now, we test if each of the five points satisfies this equation:")
    all_on_circle = True
    for i in range(len(point_list)):
        p = point_list[i]
        name = point_names[i]
        
        # Calculate the left side of the circle equation for the point p.
        lhs = (p[0] - h)**2 + (p[1] - k)**2
        
        print(f"Testing point {name}{p}:")
        print(f"  ({p[0]} - {h})^2 + ({p[1]} - {k})^2 = {lhs}")

        # Check if the result equals r_squared.
        if not math.isclose(lhs, r_squared):
            print(f"  Result {lhs} != {r_squared}. Point {name} is NOT on the circle.")
            all_on_circle = False
        else:
            print(f"  Result {lhs} == {r_squared}. Point {name} IS on the circle.")
        print()

    # Step 5: State the final conclusion based on the check.
    print("-" * 40)
    if all_on_circle:
        print("Conclusion: All five points are concyclic and can lie on a sphere.")
    else:
        print("Conclusion: The five points are not concyclic.")
        print("Since coplanar points must be concyclic to lie on a sphere, it is impossible")
        print("for all five legs to touch the surface of a sphere simultaneously.")
        print("\nTherefore, the number of such locations is 0.")

solve()