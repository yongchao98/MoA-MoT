import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface defined by three points X, Y, and Z.

    The method involves:
    1. Defining the 3D coordinates of the points based on the map.
    2. Creating two vectors on the plane (v_YX and v_YZ).
    3. Finding the normal vector to the plane via cross product.
    4. Calculating the plane's slope from the normal vector components.
    5. Determining the dip angle as the arctangent of the slope.
    """
    print("Step 1: Determine the 3D coordinates of points X, Y, and Z.")
    
    # From the map, the scale is 200 meters for 194 pixels.
    # Scale factor = 200 / 194 m/pixel
    scale = 200.0 / 194.0

    # Pixel coordinates measured from the image, with (0,0) at the top-left.
    # Y_pix = (282, 608), X_pix = (715, 523), Z_pix = (428, 290)
    # We set Y as the origin (0,0) for the map coordinates.
    # The map's y-axis points North (up), so it's inverted relative to pixel y-axis.
    x_X = (715 - 282) * scale
    y_X = (608 - 523) * scale
    x_Z = (428 - 282) * scale
    y_Z = (608 - 290) * scale

    # 3D coordinates (x, y, height) for each point.
    p_X = np.array([x_X, y_X, 120.0])
    p_Y = np.array([0.0, 0.0, 80.0])
    p_Z = np.array([x_Z, y_Z, 140.0])

    print(f"Point X: ({p_X[0]:.2f} m, {p_X[1]:.2f} m, {p_X[2]:.2f} m)")
    print(f"Point Y: ({p_Y[0]:.2f} m, {p_Y[1]:.2f} m, {p_Y[2]:.2f} m)")
    print(f"Point Z: ({p_Z[0]:.2f} m, {p_Z[1]:.2f} m, {p_Z[2]:.2f} m)\n")

    print("Step 2: Define two vectors on the plane.")
    v_YX = p_X - p_Y
    v_YZ = p_Z - p_Y
    print(f"Vector YX = X - Y = ({v_YX[0]:.2f}, {v_YX[1]:.2f}, {v_YX[2]:.2f})")
    print(f"Vector YZ = Z - Y = ({v_YZ[0]:.2f}, {v_YZ[1]:.2f}, {v_YZ[2]:.2f})\n")

    print("Step 3: Calculate the normal vector to the plane using the cross product.")
    normal_vector = np.cross(v_YX, v_YZ)
    nx, ny, nz = normal_vector
    print(f"Normal Vector n = YX x YZ = ({nx:.2f}, {ny:.2f}, {nz:.2f})\n")

    print("Step 4: Calculate the slope of the plane.")
    # Slope = sqrt(nx^2 + ny^2) / |nz|
    horizontal_component = np.sqrt(nx**2 + ny**2)
    vertical_component = np.abs(nz)
    slope = horizontal_component / vertical_component
    print(f"Slope = sqrt({nx:.2f}² + {ny:.2f}²) / |{nz:.2f}|")
    print(f"      = {horizontal_component:.2f} / {vertical_component:.2f} = {slope:.4f}\n")

    print("Step 5: Calculate the dip angle.")
    # Dip angle = arctan(slope)
    dip_rad = np.arctan(slope)
    dip_deg = np.rad2deg(dip_rad)
    print(f"Dip Angle = arctan({slope:.4f}) = {dip_deg:.2f} degrees\n")

    print("Step 6: Round the dip angle to the nearest degree.")
    rounded_dip = round(dip_deg)
    print(f"The dip of the planar surface, rounded to the nearest degree, is {rounded_dip} degrees.")
    
    return rounded_dip

if __name__ == '__main__':
    final_answer = calculate_dip()
    print(f"\n<<<9>>>")
