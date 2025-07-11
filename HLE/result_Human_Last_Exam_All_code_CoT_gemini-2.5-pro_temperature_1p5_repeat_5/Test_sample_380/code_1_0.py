import numpy as np

def explain_and_transform_coordinates():
    """
    Demonstrates the transformation from a standard coordinate system (like RTN)
    to a Field-Aligned Coordinate (FAC) system for wave analysis.
    """
    # 1. Define sample data in the original RTN (Radial-Tangential-Normal) system.
    # The local Interplanetary Magnetic Field (IMF) at L1 is not purely radial.
    # It follows the Parker Spiral. Let's assume a typical direction (in nanotesla, nT).
    # R is radial from Sun, T is tangential, N is normal to ecliptic plane.
    B0_rtn = np.array([-4.0, 2.5, 0.5])  # Background magnetic field vector

    # A fluctuation vector from an AIC wave measured at the same time.
    dB_rtn = np.array([0.1, 0.8, -0.4])  # Wave magnetic field perturbation vector

    print("--- Initial Data in Standard (RTN) Coordinates ---")
    print(f"Background Magnetic Field B0 = {B0_rtn} nT")
    print(f"Wave Perturbation dB = {dB_rtn} nT\n")
    print("The goal is to find the components of dB that are perpendicular to B0.")
    print("This requires transforming into a Field-Aligned Coordinate (FAC) system.\n")


    # 2. Create the basis vectors for the new Field-Aligned Coordinate (FAC) system.
    # The new 'z' axis is parallel to the background magnetic field B0.
    # We normalize B0 to get a unit vector.
    b_parallel_unit_vec = B0_rtn / np.linalg.norm(B0_rtn)

    # The new 'y' axis must be perpendicular to B0. We can find it using the
    # cross product with a known direction, like the original 'N' axis [0, 0, 1].
    # This choice defines the orientation of the perpendicular plane.
    b_perp1_unit_vec = np.cross(b_parallel_unit_vec, np.array([0, 0, 1]))
    b_perp1_unit_vec = b_perp1_unit_vec / np.linalg.norm(b_perp1_unit_vec)

    # The new 'x' axis completes the right-handed system (perp1, perp2, parallel).
    # Here, our perp2 vector is b_perp1 and parallel is b_parallel_unit_vec.
    # So the remaining perpendicular vector is cross(perp1, parallel).
    b_perp2_unit_vec = np.cross(b_parallel_unit_vec, b_perp1_unit_vec)

    print("--- Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    print(f"Parallel axis (z_fac):      {np.round(b_parallel_unit_vec, 3)}")
    print(f"Perpendicular axis 1 (x_fac): {np.round(b_perp2_unit_vec, 3)}")
    print(f"Perpendicular axis 2 (y_fac): {np.round(b_perp1_unit_vec, 3)}\n")

    # 3. Project the wave vector (dB_rtn) onto the new FAC basis vectors
    #    to find its components in the new system.
    db_comp_parallel = np.dot(dB_rtn, b_parallel_unit_vec)
    db_comp_perp1 = np.dot(dB_rtn, b_perp2_unit_vec)
    db_comp_perp2 = np.dot(dB_rtn, b_perp1_unit_vec)

    print("--- Wave Components in the Field-Aligned System ---")
    print("These are the physically meaningful components for helicity calculation.")
    print(f"Component of dB parallel to B0:       {db_comp_parallel:.4f} nT")
    print(f"Component of dB perpendicular to B0 (1): {db_comp_perp1:.4f} nT")
    print(f"Component of dB perpendicular to B0 (2): {db_comp_perp2:.4f} nT\n")

    print("CONCLUSION:")
    print("The normalized magnetic helicity would be calculated using a time series of the two")
    print(f"perpendicular components: {db_comp_perp1:.4f} and {db_comp_perp2:.4f}, not the original T and N components.")

if __name__ == '__main__':
    explain_and_transform_coordinates()