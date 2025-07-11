import numpy as np

def calculate_field_aligned_components():
    """
    This script demonstrates the correct way to find the magnetic field components
    perpendicular to the local magnetic field for helicity calculations.
    """
    # Let's assume we are in a standard coordinate system like GSE (Geocentric Solar Ecliptic)
    # where X is roughly radial from the Sun, Y is duskward, and Z is ecliptic north.
    
    # 1. Define the average local background magnetic field (B0) at L1.
    # Note that it is NOT purely radial (i.e., not just in X). Due to the Parker Spiral,
    # it has a significant Y component. Units are nanoTesla (nT).
    # Bx = -4.5 nT (points toward Sun), By = 3.5 nT, Bz = -1.0 nT
    B0 = np.array([-4.5, 3.5, -1.0])

    # 2. Define a sample magnetic field fluctuation vector (b_wave) at a single point in time.
    # This represents the perturbation caused by the AIC wave.
    b_wave = np.array([0.1, 0.4, -0.2])

    print("--- Initial Data in GSE-like Coordinates ---")
    print(f"Average Local Magnetic Field (B0): {B0} nT")
    print(f"Wave Fluctuation (b_wave): {b_wave} nT")
    print(f"Component perpendicular to radial X-axis (by): {b_wave[1]}")
    print(f"Component perpendicular to radial X-axis (bz): {b_wave[2]}\n")

    # 3. Create the Field-Aligned Coordinate (FAC) system based on B0.
    
    # The parallel unit vector (ez_fac) is in the direction of B0.
    ez_fac = B0 / np.linalg.norm(B0)

    # To create the perpendicular axes, we take the cross product with a fixed vector.
    # We use the Z-axis of the original system, k_hat = [0, 0, 1].
    # This choice is arbitrary but standard. If B0 is parallel to k_hat, we would choose another axis.
    k_hat = np.array([0., 0., 1.])
    
    # One perpendicular axis (ey_fac) is ez_fac x k_hat, then normalized.
    ey_fac = np.cross(ez_fac, k_hat)
    ey_fac = ey_fac / np.linalg.norm(ey_fac)

    # The other perpendicular axis (ex_fac) completes the right-handed system.
    ex_fac = np.cross(ey_fac, ez_fac)

    print("--- Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    print(f"Parallel axis (ez_fac, along B0): {np.round(ez_fac, 3)}")
    print(f"Perp. axis 1 (ex_fac):             {np.round(ex_fac, 3)}")
    print(f"Perp. axis 2 (ey_fac):             {np.round(ey_fac, 3)}\n")

    # 4. Project the wave vector onto the new perpendicular FAC axes.
    # This is done using the dot product.
    b_perp1 = np.dot(b_wave, ex_fac)
    b_perp2 = np.dot(b_wave, ey_fac)
    b_parallel = np.dot(b_wave, ez_fac)

    print("--- Wave Components in Field-Aligned Coordinates ---")
    print("These are the physically correct components for helicity calculation.")
    print(f"Component 1 perpendicular to B0 (b_perp1): {b_perp1:.4f} nT")
    print(f"Component 2 perpendicular to B0 (b_perp2): {b_perp2:.4f} nT")
    print(f"Component parallel to B0 (b_parallel):     {b_parallel:.4f} nT\n")
    
    print("--- Justification ---")
    print("As you suspected, the local magnetic field at L1 is generally NOT radial due to the Parker Spiral.")
    print("Therefore, for a physically meaningful helicity calculation of waves that propagate along the magnetic field,")
    print("one MUST transform the data into a Field-Aligned Coordinate (FAC) system.")
    print("The normalized magnetic helicity should be calculated using the Fourier transforms of the two perpendicular")
    print(f"components found above (b_perp1 = {b_perp1:.4f} and b_perp2 = {b_perp2:.4f}), NOT the original by and bz components.")

if __name__ == '__main__':
    calculate_field_aligned_components()
