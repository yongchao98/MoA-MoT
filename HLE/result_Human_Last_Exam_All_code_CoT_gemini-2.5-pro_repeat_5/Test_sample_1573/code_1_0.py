import math

def main():
    """
    Checks if the five points of the chair legs are concyclic.
    """
    # The five points at the base of the chair legs
    points = [(0, 0), (2, 0), (2, 2), (0, 2), (1, 4)]
    
    # The first four points form a square. The center of their circumcircle is the
    # center of the square.
    center_x = (points[0][0] + points[2][0]) / 2.0
    center_y = (points[0][1] + points[2][1]) / 2.0
    center = (center_x, center_y)

    # The radius-squared is the squared distance from the center to any of the first four points.
    # We use the first point (0,0).
    radius_sq = (points[0][0] - center_x)**2 + (points[0][1] - center_y)**2
    
    print("Step 1: Determine the circle defined by the first four legs.")
    print(f"The first four legs are at {points[:4]}.")
    print(f"These points form a square. The center of the circumcircle is {center}.")
    print(f"The equation of a circle is (x - h)^2 + (y - k)^2 = r^2.")
    print(f"The equation for the circle through the first four points is (x - {int(center_x)})^2 + (y - {int(center_y)})^2 = {int(radius_sq)}.\n")

    # Now, check if the fifth point lies on this circle.
    fifth_point = points[4]
    lhs = (fifth_point[0] - center_x)**2 + (fifth_point[1] - center_y)**2
    
    print("Step 2: Check if the fifth leg's point (1, 4) is on this circle.")
    print(f"Plugging the fifth point {fifth_point} into the left side of the equation:")
    print(f"({fifth_point[0]} - {int(center_x)})^2 + ({fifth_point[1]} - {int(center_y)})^2 = {int(lhs)}")
    
    print("\nStep 3: Compare the result with the required radius-squared.")
    print(f"The result is {int(lhs)}, but the radius-squared of the circle is {int(radius_sq)}.")

    if math.isclose(lhs, radius_sq):
        print("\nConclusion: The five points are concyclic.")
        print("This means there is an infinite number of places the chair could fit.")
    else:
        print("\nConclusion: The five points are NOT concyclic.")
        print("Therefore, it is impossible for all five legs to touch the surface of a sphere simultaneously.")
        print("The set of such locations is empty, and its cardinality is 0.")

if __name__ == "__main__":
    main()
