import numpy as np

def calculate_dip():
    """
    Calculates the dip of a planar surface given three points on a map.
    """
    # Step 1: Define points and their properties based on the map.
    # Pixel coordinates are measured from the image (top-left origin).
    px_coords = {
        'X': (450, 408),
        'Y': (195, 469),
        'Z': (324, 190)
    }
    # Heights are given in meters.
    heights = {
        'X': 120,
        'Y': 80,
        'Z': 140
    }
    # Map scale: 100 meters corresponds to 80 pixels.
    scale_m_per_px = 100.0 / 80.0
    # Image height in pixels to invert y-axis for standard Cartesian coordinates.
    image_height_px = 774

    print("Step 1: Convert map locations to 3D coordinates (x, y, z).")
    
    # Create 3D points in meters.
    # x_map = px * scale
    # y_map = (image_height - py) * scale
    p_x = np.array([
        px_coords['X'][0] * scale_m_per_px,
        (image_height_px - px_coords['X'][1]) * scale_m_per_px,
        heights['X']
    ])
    p_y = np.array([
        px_coords['Y'][0] * scale_m_per_px,
        (image_height_px - px_coords['Y'][1]) * scale_m_per_px,
        heights['Y']
    ])
    p_z = np.array([
        px_coords['Z'][0] * scale_m_per_px,
        (image_height_px - px_coords['Z'][1]) * scale_m_per_px,
        heights['Z']
    ])

    print(f"  - Point X: ({p_x[0]:.2f} m, {p_x[1]:.2f} m, {p_x[2]:.2f} m)")
    print(f"  - Point Y: ({p_y[0]:.2f} m, {p_y[1]:.2f} m, {p_y[2]:.2f} m)")
    print(f"  - Point Z: ({p_z[0]:.2f} m, {p_z[1]:.2f} m, {p_z[2]:.2f} m)\n")

    # Step 2: Create two vectors on the plane.
    vec_yx = p_x - p_y
    vec_yz = p_z - p_y
    print("Step 2: Define two vectors on the plane (e.g., YX and YZ).")
    print(f"  - Vector YX (X - Y): [{vec_yx[0]:.2f}, {vec_yx[1]:.2f}, {vec_yx[2]:.2f}]")
    print(f"  - Vector YZ (Z - Y): [{vec_yz[0]:.2f}, {vec_yz[1]:.2f}, {vec_yz[2]:.2f}]\n")

    # Step 3: Calculate the normal vector to the plane.
    normal_vector = np.cross(vec_yx, vec_yz)
    nx, ny, nz = normal_vector
    print("Step 3: Calculate the normal vector N by taking the cross product of YX and YZ.")
    print(f"  - Normal Vector N = [{nx:.2f}, {ny:.2f}, {nz:.2f}]\n")

    # Step 4: Calculate the dip angle.
    horizontal_component_mag = np.sqrt(nx**2 + ny**2)
    vertical_component_mag = np.abs(nz)
    tan_dip = horizontal_component_mag / vertical_component_mag
    dip_rad = np.arctan(tan_dip)
    dip_deg = np.degrees(dip_rad)
    
    print("Step 4: Calculate the dip angle (δ) using the normal vector components.")
    print("  - The formula for the tangent of the dip is: tan(δ) = sqrt(Nx² + Ny²) / |Nz|")
    print(f"  - tan(δ) = sqrt({nx:.2f}² + {ny:.2f}²) / |{nz:.2f}|")
    print(f"  - tan(δ) = {horizontal_component_mag:.2f} / {vertical_component_mag:.2f} = {tan_dip:.4f}")
    print(f"  - δ = arctan({tan_dip:.4f}) = {dip_deg:.2f} degrees\n")
    
    # Step 5: Round the result.
    rounded_dip = round(dip_deg)
    print("Step 5: Round the dip angle to the nearest degree.")
    print(f"  - The final dip is {rounded_dip} degrees.")
    
    return rounded_dip

# Execute the calculation and print the final answer in the required format.
final_answer = calculate_dip()
# The final answer needs to be enclosed in <<<>>>
# print(f"\n<<< {final_answer} >>>")
# Since the required output format doesn't want the function's return value to be printed in the prompt, let's just make the function do all the printing and have a final statement
print(f"\nThe dip of the planar surface, rounded to the nearest degree, is {final_answer} degrees.")
print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    # The script is self-contained. Calling calculate_dip() does everything.
    pass
