import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface from three points on a map.
    """
    # Step 1: Define coordinates and scale from the map analysis.
    # The North arrow is vertical, so the pixel y-axis is aligned with the North-South direction.
    # We assume the origin (0,0) of the pixel coordinate system is at the top-left of the image.
    
    # Pixel coordinates (x, y)
    X_pix = np.array([450, 345])
    Y_pix = np.array([145, 425])
    Z_pix = np.array([330, 140])

    # Heights in meters
    h_X = 120
    h_Y = 80
    h_Z = 140

    # The scale bar shows that 200 meters correspond to approximately 158 pixels.
    m_per_px = 200.0 / 158.0
    
    print("Step 1: Establishing a real-world coordinate system (x=East, y=North, z=Up).")
    print(f"Point Y is at height {h_Y} m.")
    print(f"Point X is at height {h_X} m.")
    print(f"Point Z is at height {h_Z} m.")
    print(f"Map scale: {m_per_px:.4f} meters/pixel.\n")

    # Step 2: Set Y as the origin of the (x, y) plane and calculate coordinates for X and Z.
    # The positive y-direction is North, which corresponds to a decreasing pixel y-coordinate.
    x_X_m = (X_pix[0] - Y_pix[0]) * m_per_px
    y_X_m = (Y_pix[1] - X_pix[1]) * m_per_px  # Invert y-axis for North
    h_X_rel = h_X - h_Y

    x_Z_m = (Z_pix[0] - Y_pix[0]) * m_per_px
    y_Z_m = (Y_pix[1] - Z_pix[1]) * m_per_px  # Invert y-axis for North
    h_Z_rel = h_Z - h_Y
    
    print("Step 2: Defining the plane with equation z = Ax + By + C.")
    print("By placing point Y at the origin, the equation for relative height is z' = Ax + By.")
    print("We need to solve the following system of equations:")
    print(f"{x_X_m:.1f}*A + {y_X_m:.1f}*B = {h_X_rel}")
    print(f"{x_Z_m:.1f}*A + {y_Z_m:.1f}*B = {h_Z_rel}\n")

    # Step 3: Solve the system of linear equations for A and B.
    # [x_X, y_X] [A] = [h_X_rel]
    # [x_Z, y_Z] [B]   [h_Z_rel]
    coeffs = np.array([
        [x_X_m, y_X_m],
        [x_Z_m, y_Z_m]
    ])
    vals = np.array([h_X_rel, h_Z_rel])
    
    try:
        A, B = np.linalg.solve(coeffs, vals)
        print(f"Step 3: Solved for the gradient components.")
        print(f"A = {A:.4f}")
        print(f"B = {B:.4f}\n")
    except np.linalg.LinAlgError:
        print("The points are collinear and do not define a unique plane.")
        return

    # Step 4: Calculate the dip.
    # The tangent of the dip angle is the magnitude of the gradient.
    gradient_magnitude = np.sqrt(A**2 + B**2)
    dip_radians = np.arctan(gradient_magnitude)
    dip_degrees = np.degrees(dip_radians)

    print("Step 4: Calculating the dip angle.")
    print("The dip angle is calculated as: arctan(sqrt(A^2 + B^2))")
    print(f"dip = arctan(sqrt({A:.4f}^2 + {B:.4f}^2))")
    print(f"dip = arctan({gradient_magnitude:.4f})")
    print(f"The calculated dip is {dip_degrees:.2f} degrees.\n")

    # Step 5: Round the result to the nearest degree.
    rounded_dip = int(round(dip_degrees))
    print(f"Step 5: Rounding to the nearest degree.")
    print(f"The dip of the planar surface is {rounded_dip} degrees.")

calculate_dip()