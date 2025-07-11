import numpy as np

def calculate_field_aligned_components():
    """
    Demonstrates the transformation of magnetic field wave components into a
    Field-Aligned Coordinate (FAC) system for accurate helicity calculation.
    """
    # Step 1: Define a realistic background magnetic field at L1 in GSE coordinates (nT).
    # The field is not radial due to the Parker Spiral. Bx is negative (points from Earth to Sun),
    # and By is significant. Bz is typically small.
    # B0_gse = [Bx, By, Bz]
    B0_gse = np.array([-4.5, 4.0, 0.5])

    # Step 2: Define a sample magnetic field fluctuation (the AIC wave) in GSE coordinates (nT).
    delta_B_gse = np.array([0.3, -0.2, 0.6])
    
    print("--- Input Data (in GSE coordinates) ---")
    print(f"Background Magnetic Field (B0): {B0_gse} nT")
    print(f"Wave Fluctuation (δB): {delta_B_gse} nT\n")

    # --- Incorrect Approach ---
    # Using components perpendicular to the radial (X) direction directly.
    delta_B_y_gse = delta_B_gse[1]
    delta_B_z_gse = delta_B_gse[2]
    print("--- Simplified (Incorrect) Approach ---")
    print("Using components perpendicular to the radial (X) direction:")
    print(f"δB_perp1 (δBy_gse) = {delta_B_y_gse:.3f} nT")
    print(f"δB_perp2 (δBz_gse) = {delta_B_z_gse:.3f} nT\n")

    # --- Correct Approach: Transform to Field-Aligned Coordinates (FAC) ---
    print("--- Physically Correct Approach ---")
    print("Transforming to a Field-Aligned Coordinate (FAC) system:\n")
    
    # Step 3: Create the FAC system basis vectors.
    # The parallel unit vector is along the background magnetic field B0.
    b_parallel = B0_gse / np.linalg.norm(B0_gse)

    # To create the perpendicular vectors, we use cross products.
    # We need a reference vector that is not parallel to B0. The GSE Y-axis [0,1,0] is a good choice.
    ref_vector = np.array([0.0, 1.0, 0.0])
    
    # The first perpendicular axis (b_perp2) is B0 x Y_gse, normalized.
    # This creates a vector perpendicular to the plane containing the magnetic field and the Y-axis.
    b_perp2 = np.cross(b_parallel, ref_vector)
    b_perp2 = b_perp2 / np.linalg.norm(b_perp2)

    # The second perpendicular axis (b_perp1) completes the right-hand system.
    b_perp1 = np.cross(b_perp2, b_parallel)

    print("FAC Unit Vectors:")
    print(f"b_parallel (along B0)  = {np.round(b_parallel, 3)}")
    print(f"b_perp1                = {np.round(b_perp1, 3)}")
    print(f"b_perp2                = {np.round(b_perp2, 3)}\n")

    # Step 4: Project the wave fluctuation (delta_B_gse) onto the new perpendicular basis vectors.
    delta_B_perp1_val = np.dot(delta_B_gse, b_perp1)
    delta_B_perp2_val = np.dot(delta_B_gse, b_perp2)

    print("Correct Perpendicular Wave Components in FAC system:")
    print(f"δB_perp1 = {delta_B_perp1_val:.3f} nT")
    print(f"δB_perp2 = {delta_B_perp2_val:.3f} nT\n")

    # Step 5: Show the final conceptual equation for helicity.
    # The normalized magnetic helicity, σ_m, is calculated from the cross-spectrum of these
    # two perpendicular components over a time series. The equation looks conceptually like:
    # σ_m = Im(<δB_perp1(f) * δB_perp2(f)*>) / (<|δB_perp1(f)|²> + <|δB_perp2(f)|²>)
    # For this single data point, we can show where the numbers would be plugged in.
    
    print("--- Final Equation for Helicity ---")
    print("The normalized magnetic helicity requires a time series, but its calculation relies on the correct components.")
    print("A key term in its calculation, which determines the polarization, is built from these components.")
    print("Using the calculated values, the term for a single time point looks like:")
    print(f"Polarization Term ∝ (δB_perp1 * H(δB_perp2) - δB_perp2 * H(δB_perp1))")
    print(f"Which uses the numbers: ({delta_B_perp1_val:.3f} * H({delta_B_perp2_val:.3f})) - ({delta_B_perp2_val:.3f} * H({delta_B_perp1_val:.3f}))")
    print("(where H is the Hilbert Transform, representing a 90-degree phase shift)")

if __name__ == '__main__':
    calculate_field_aligned_components()
