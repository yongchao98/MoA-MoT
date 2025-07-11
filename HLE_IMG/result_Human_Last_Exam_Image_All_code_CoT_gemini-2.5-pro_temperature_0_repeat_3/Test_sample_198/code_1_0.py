import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface from three points with known coordinates and heights.
    """
    # Step 1: Define the coordinates of the three points.
    # The (x, y) coordinates are estimated from the map image, where the scale is 1 pixel = 1 meter.
    # The heights (z) are given in the problem description.
    # Point coordinates (x, y, z) in meters:
    X = np.array([553, 433, 120])
    Y = np.array([203, 283, 80])
    Z = np.array([403, 583, 140])

    print("Step 1: Define the 3D coordinates of points X, Y, and Z.")
    print(f"Point X: (x={X[0]} m, y={X[1]} m, height={X[2]} m)")
    print(f"Point Y: (x={Y[0]} m, y={Y[1]} m, height={Y[2]} m)")
    print(f"Point Z: (x={Z[0]} m, y={Z[1]} m, height={Z[2]} m)")
    print("-" * 30)

    # Step 2: Create two vectors lying on the planar surface.
    # We'll use vectors originating from point Y for this calculation.
    vec_YX = X - Y
    vec_YZ = Z - Y

    print("Step 2: Create two vectors on the plane (e.g., from Y to X and Y to Z).")
    print(f"Vector YX = X - Y = <{vec_YX[0]}, {vec_YX[1]}, {vec_YX[2]}>")
    print(f"Vector YZ = Z - Y = <{vec_YZ[0]}, {vec_YZ[1]}, {vec_YZ[2]}>")
    print("-" * 30)

    # Step 3: Calculate the normal vector to the plane using the cross product.
    normal_vector = np.cross(vec_YX, vec_YZ)

    print("Step 3: Calculate the normal vector (n) to the plane using the cross product.")
    print(f"n = YX x YZ = <{normal_vector[0]}, {normal_vector[1]}, {normal_vector[2]}>")
    print("-" * 30)

    # Step 4: Calculate the dip angle.
    # The dip is the angle between the plane's normal vector and the vertical axis (0,0,1).
    vertical_vector = np.array([0, 0, 1])
    
    # Using the dot product formula: angle = arccos((v1 . v2) / (|v1| * |v2|))
    dot_product = np.dot(normal_vector, vertical_vector)
    mag_normal = np.linalg.norm(normal_vector)
    mag_vertical = np.linalg.norm(vertical_vector) # This is always 1

    # We use the absolute value of the dot product to ensure we get the acute angle.
    cos_angle = abs(dot_product) / (mag_normal * mag_vertical)
    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)

    print("Step 4: Calculate the dip angle (δ).")
    print("This is the angle between the normal vector and the vertical vector k = <0, 0, 1>.")
    print(f"The equation is: δ = arccos( |n · k| / (||n|| * ||k||) )")
    print(f"δ = arccos( |{dot_product:.0f}| / ({mag_normal:.2f} * {mag_vertical:.1f}) )")
    print(f"δ = arccos({cos_angle:.4f})")
    print(f"δ = {angle_rad:.4f} radians")
    print(f"δ = {angle_deg:.2f} degrees")
    print("-" * 30)

    # Step 5: Round the answer to the nearest degree.
    rounded_dip = round(angle_deg)

    print("Step 5: Round the result to the nearest degree.")
    print(f"The final calculated dip is {rounded_dip} degrees.")

if __name__ == '__main__':
    calculate_dip()