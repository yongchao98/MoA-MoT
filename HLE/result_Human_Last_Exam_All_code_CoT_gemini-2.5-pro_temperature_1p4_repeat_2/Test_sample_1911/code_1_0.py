import numpy as np

def check_for_corners(points, tolerance_deg=30):
    """
    Checks a list of 2D points for corners. A corner is detected if the angle
    formed by a point and its two neighbors deviates significantly from 180 degrees.
    """
    corner_points = []
    # We need at least 3 points to check for a corner
    if len(points) < 3:
        return corner_points

    # Sort points by x-coordinate to ensure neighbors are adjacent in the list
    sorted_points = sorted(points, key=lambda p: p[0])

    for i in range(1, len(sorted_points) - 1):
        p_prev = np.array(sorted_points[i-1])
        p_curr = np.array(sorted_points[i])
        p_next = np.array(sorted_points[i+1])

        # Create vectors from the current point to its neighbors
        v1 = p_prev - p_curr
        v2 = p_next - p_curr

        # Normalize the vectors to unit vectors
        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)

        # Calculate the dot product, which is cos(theta) for unit vectors
        dot_product = np.dot(v1_u, v2_u)
        
        # Clip the dot product to avoid floating point errors with arccos
        dot_product = np.clip(dot_product, -1.0, 1.0)
        
        # Calculate angle in degrees
        angle_deg = np.rad2deg(np.arccos(dot_product))

        # A smooth curve will have an angle close to 180 degrees.
        # A corner will have an angle significantly different from 180.
        if abs(angle_deg - 180) > tolerance_deg:
            corner_points.append(sorted_points[i])
            
    return corner_points

def main():
    # 1. Define the set L by sampling points
    x_vals = np.linspace(-5, 5, 101)
    L = [(x, abs(x)) for x in x_vals]

    # 2. Case 1: Remove the point z = (0,0)
    z1 = (0.0, 0.0)
    L_without_z1 = [p for p in L if p != z1]
    
    # Check for corners in L_without_z1.
    # We need to split the remaining set into its two branches to check corners correctly.
    left_branch = [p for p in L_without_z1 if p[0] < 0]
    right_branch = [p for p in L_without_z1 if p[0] > 0]
    
    corners1_left = check_for_corners(left_branch)
    corners1_right = check_for_corners(right_branch)
    total_corners1 = corners1_left + corners1_right
    
    print(f"Testing statement E computationally.")
    print("-" * 30)
    print(f"Case 1: Removing z = {z1}")
    if not total_corners1:
        print("Result: No corners found. The remaining set is smooth.")
    else:
        print(f"Result: Found corner(s) at {total_corners1}")

    # 3. Case 2: Remove a point z != (0,0), for example z = (2.0, 2.0)
    # Find the point in L closest to (2.0, 2.0) to remove
    z2_x = 2.0
    z2 = (z2_x, abs(z2_x))
    L_without_z2 = [p for p in L if p != z2]
    
    corners2 = check_for_corners(L_without_z2)
    
    print("-" * 30)
    print(f"Case 2: Removing z = {z2}")
    if not corners2:
        print("Result: No corners found.")
    else:
        # We expect to find the corner at (0,0)
        detected_corner_x = corners2[0][0]
        detected_corner_y = corners2[0][1]
        print(f"Result: Found a corner at ({detected_corner_x:.1f}, {detected_corner_y:.1f})")

main()
