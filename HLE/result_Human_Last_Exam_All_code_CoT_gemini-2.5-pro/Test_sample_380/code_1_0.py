import numpy as np

def solve_task():
    """
    Demonstrates the rotation of magnetic field data into a Field-Aligned
    Coordinate (FAC) system to justify the calculation of helicity.
    """
    print("This script demonstrates why 'Y and Z components' are used to calculate helicity for AIC waves.")
    print("It's not because the local field is radial, but because the data is rotated into a Field-Aligned Coordinate (FAC) system first.\n")

    # 1. Define the average background magnetic field (B0) in a GSE-like system (X=radial).
    # This vector represents a typical Parker Spiral at 1 AU, which is NOT radial.
    # The X and Y components are strong, representing the spiral.
    B0_gse = np.array([-4.5, -4.0, 0.5])
    print(f"Step 1: Define the local average magnetic field B0 in the original (GSE) frame.")
    print(f"B0_gse = [{B0_gse[0]:.2f}, {B0_gse[1]:.2f}, {B0_gse[2]:.2f}] (in nT)")
    print("Note: This vector is not aligned with the X-axis (the radial direction).\n")

    # 2. Define the new basis vectors for the FAC system.
    # The new x-axis (x_fac) is parallel to B0.
    x_fac = B0_gse / np.linalg.norm(B0_gse)

    # The new y-axis (y_fac) is defined to be perpendicular to B0 and the GSE Z-axis.
    # This keeps the new y-axis in the ecliptic plane, a common convention.
    z_gse_axis = np.array([0., 0., 1.])
    y_fac = np.cross(z_gse_axis, x_fac)
    y_fac = y_fac / np.linalg.norm(y_fac)

    # The new z-axis (z_fac) completes the right-handed system.
    z_fac = np.cross(x_fac, y_fac)

    print("Step 2: Construct the new FAC system basis vectors based on B0.")
    print(f"x_fac (parallel to B0)  = [{x_fac[0]:.3f}, {x_fac[1]:.3f}, {x_fac[2]:.3f}]")
    print(f"y_fac (perpendicular 1) = [{y_fac[0]:.3f}, {y_fac[1]:.3f}, {y_fac[2]:.3f}]")
    print(f"z_fac (perpendicular 2) = [{z_fac[0]:.3f}, {z_fac[1]:.3f}, {z_fac[2]:.3f}]\n")

    # 3. Define a sample magnetic field fluctuation vector (deltaB) in the original GSE frame.
    # For an AIC wave, this fluctuation should be nearly perpendicular to B0.
    # We create a fluctuation that lives mostly in the new y_fac-z_fac plane.
    # deltaB_gse = 0.2 * y_fac + 0.3 * z_fac
    deltaB_gse = 0.2 * y_fac + 0.3 * z_fac
    print("Step 3: Define a sample magnetic field fluctuation vector, deltaB, in the original GSE frame.")
    print(f"deltaB_gse = [{deltaB_gse[0]:.3f}, {deltaB_gse[1]:.3f}, {deltaB_gse[2]:.3f}] (in nT)\n")
    
    # Check that deltaB_gse is indeed perpendicular to B0_gse
    dot_product_check = np.dot(deltaB_gse, B0_gse)
    print(f"Check: The dot product of deltaB_gse and B0_gse is {dot_product_check:.4f}, which is close to zero, confirming they are perpendicular.\n")


    # 4. Rotate the deltaB vector into the new FAC system.
    # The components of the new vector are the dot products of the original vector with the new basis vectors.
    deltaB_fac_x = np.dot(deltaB_gse, x_fac)
    deltaB_fac_y = np.dot(deltaB_gse, y_fac)
    deltaB_fac_z = np.dot(deltaB_gse, z_fac)

    print("Step 4: Rotate deltaB into the FAC system to find its parallel and perpendicular components.")
    print("The final equation for the rotated vector is:")
    print(f"deltaB_fac = [{deltaB_fac_x:.3f}, {deltaB_fac_y:.3f}, {deltaB_fac_z:.3f}]")
    print(f"Parallel component (x'): {deltaB_fac_x:.3f}")
    print(f"Perpendicular component (y'): {deltaB_fac_y:.3f}")
    print(f"Perpendicular component (z'): {deltaB_fac_z:.3f}\n")

    print("Conclusion:")
    print("The normalized magnetic helicity would now be calculated using the time series of these perpendicular components (y' and z').")
    print("This shows the use of 'Y and Z components' is justified because the analysis is performed in a coordinate system aligned with the local magnetic field.")

solve_task()