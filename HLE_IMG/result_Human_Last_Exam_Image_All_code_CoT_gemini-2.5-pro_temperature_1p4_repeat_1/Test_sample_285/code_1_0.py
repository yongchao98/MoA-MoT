import math

def get_circumcircle(p1, p2, p3):
    """
    Calculates the circumcenter and the square of the circumradius of a triangle.
    p1, p2, p3 are tuples (x, y).
    Returns a tuple: ((center_x, center_y), radius_squared)
    """
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    
    if D == 0:
        # Points are collinear, circumcircle is undefined (infinite radius)
        return (None, float('inf'))

    # Calculate center_x
    ux = ((p1[0]**2 + p1[1]**2) * (p2[1] - p3[1]) + 
          (p2[0]**2 + p2[1]**2) * (p3[1] - p1[1]) + 
          (p3[0]**2 + p3[1]**2) * (p1[1] - p2[1])) / D
          
    # Calculate center_y
    uy = ((p1[0]**2 + p1[1]**2) * (p3[0] - p2[0]) + 
          (p2[0]**2 + p2[1]**2) * (p1[0] - p3[0]) + 
          (p3[0]**2 + p3[1]**2) * (p2[0] - p1[0])) / D
          
    center = (ux, uy)
    
    # Calculate radius squared
    radius_sq = (p1[0] - ux)**2 + (p1[1] - uy)**2
    
    return (center, radius_sq)

def is_inside_circle(point, circle_center, circle_radius_sq):
    """
    Checks if a point is inside a circle.
    """
    dist_sq = (point[0] - circle_center[0])**2 + (point[1] - circle_center[1])**2
    return dist_sq < circle_radius_sq

# Step 1: Assign approximate coordinates based on diagram D
points_d = {
    'Outer Top': (0, 10),
    'Outer Top-Right': (4, 3),
    'Inner Center': (0, 0),
    'Inner Right': (3, 2)
}

# Step 2: Define the triangle to check
p1 = points_d['Outer Top']
p2 = points_d['Outer Top-Right']
p3 = points_d['Inner Center']
triangle_to_check = (p1, p2, p3)
print(f"Checking the triangle with vertices: P1={p1}, P2={p2}, P3={p3}")

# Step 3: Define the point to test
point_to_test = points_d['Inner Right']
print(f"Testing if point P_test={point_to_test} is inside its circumcircle.")
print("-" * 20)

# Step 4: Calculate the circumcircle
circumcenter, radius_sq = get_circumcircle(p1, p2, p3)
print(f"Calculated Circumcenter: ({circumcenter[0]:.4f}, {circumcenter[1]:.4f})")
print(f"Calculated Circumradius squared: {radius_sq:.4f}")

# Step 5: Calculate the distance from the test point to the circumcenter
dist_sq_from_center = (point_to_test[0] - circumcenter[0])**2 + (point_to_test[1] - circumcenter[1])**2
print(f"Squared distance of P_test from circumcenter: {dist_sq_from_center:.4f}")
print("-" * 20)

# Step 6: Compare distance and radius
is_inside = dist_sq_from_center < radius_sq
print(f"Is the squared distance less than the squared radius? {is_inside}")

if is_inside:
    print(f"\nConclusion: Point {point_to_test} is inside the circumcircle of triangle {triangle_to_check}.")
    print("Therefore, triangulation D violates the empty circle property.")
else:
    print("\nConclusion: The chosen point is not inside the circumcircle of the chosen triangle.")
