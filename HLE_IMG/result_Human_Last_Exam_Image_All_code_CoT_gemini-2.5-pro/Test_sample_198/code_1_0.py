import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define Input Data from the map
    # Heights of the points in meters
    h_X = 120
    h_Y = 80
    h_Z = 140

    # Pixel coordinates (x, y) measured from the top-left corner of the image
    pix_X = np.array([568, 417])
    pix_Y = np.array([208, 497])
    pix_Z = np.array([369, 238])

    # Scale bar information: 200 meters corresponds to 145 pixels (795 - 650).
    scale_factor = 200.0 / 145.0

    # Step 2: Convert pixel coordinates to a 3D map coordinate system
    # We set point Y as the origin (0, 0) for the x-y plane.
    # The map's North is up, so we invert the y-axis pixel direction.
    map_X_xy = (pix_X - pix_Y) * scale_factor
    map_X_xy[1] *= -1  # Invert y-axis to point North
    
    map_Y_xy = np.array([0.0, 0.0]) # Y is the origin
    
    map_Z_xy = (pix_Z - pix_Y) * scale_factor
    map_Z_xy[1] *= -1 # Invert y-axis to point North

    # Create 3D points P = (x, y, z)
    P_X = np.append(map_X_xy, h_X)
    P_Y = np.append(map_Y_xy, h_Y)
    P_Z = np.append(map_Z_xy, h_Z)
    
    # Step 3: Define two vectors on the plane
    # Vectors start from point Y to avoid dealing with offsets.
    vec_YX = P_X - P_Y
    vec_YZ = P_Z - P_Y

    # Step 4: Find the normal vector to the plane via cross product
    normal_vector = np.cross(vec_YX, vec_YZ)
    Nx, Ny, Nz = normal_vector

    # Step 5: Calculate the dip from the normal vector components
    # The slope is the magnitude of the horizontal gradient.
    # slope = sqrt(Nx^2 + Ny^2) / |Nz|
    slope = np.sqrt(Nx**2 + Ny**2) / np.abs(Nz)

    # The dip angle is the arctangent of the slope.
    dip_degrees = np.degrees(np.arctan(slope))
    
    # Round to the nearest degree for the final answer.
    dip_rounded = round(dip_degrees)

    # Print the calculation process
    print("This script calculates the dip of the planar surface.\n")
    print(f"Point Information:")
    print(f"  X: height={h_X}m, Y: height={h_Y}m, Z: height={h_Z}m\n")
    print(f"Normal vector (N) components derived from the points: ")
    print(f"  Nx = {Nx:.2f}, Ny = {Ny:.2f}, Nz = {Nz:.2f}\n")
    print("The dip angle is calculated using the formula: Dip = arctan(sqrt(Nx^2 + Ny^2) / |Nz|)\n")
    
    # Print the final equation with numbers
    print("Final Calculation:")
    print(f"Dip = arctan(sqrt({Nx:.2f}^2 + {Ny:.2f}^2) / |{Nz:.2f}|)")
    print(f"Dip = arctan({np.sqrt(Nx**2 + Ny**2):.2f} / {np.abs(Nz):.2f})")
    print(f"Dip = arctan({slope:.4f})")
    print(f"Dip = {dip_degrees:.2f} degrees")
    print("\n----------------------------------------------------")
    print(f"The calculated dip of the planar surface is {dip_rounded} degrees.")
    print("----------------------------------------------------\n")
    return dip_rounded

# Run the calculation and store the answer
final_answer = calculate_dip()
print(f'<<<{final_answer}>>>')