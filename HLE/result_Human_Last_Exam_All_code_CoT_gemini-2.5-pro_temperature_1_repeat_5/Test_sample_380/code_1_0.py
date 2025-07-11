import numpy as np

def calculate_field_aligned_coordinates():
    """
    Demonstrates the calculation of a Field-Aligned Coordinate (FAC) system
    for a non-radial mean magnetic field, justifying the correct method for
    calculating magnetic helicity.
    """
    # Step 1: Define a realistic mean magnetic field vector at L1 (in nT, GSE coordinates).
    # The solar wind flows roughly along the -X GSE direction.
    # Due to the Parker Spiral, the magnetic field has a significant Y component.
    # A 45-degree spiral angle means |Bx| is similar to |By|.
    # B_mean = [Bx, By, Bz]
    B_mean = np.array([-4.0, -4.0, 1.0])

    # Step 2: Calculate the unit vector parallel to the local mean magnetic field.
    # This defines the direction of propagation for field-aligned waves like AICs.
    norm_B_mean = np.linalg.norm(B_mean)
    e_parallel = B_mean / norm_B_mean

    # Step 3: Construct the perpendicular plane using cross products.
    # We use a fixed reference vector, typically the radial direction (GSE X-axis),
    # to unambiguously define the new basis.
    # Reference vector (GSE X-axis, pointing from Earth towards the Sun)
    ref_vec_radial = np.array([1.0, 0.0, 0.0])

    # Step 4: Calculate the two perpendicular basis vectors.
    # e_perp2 is perpendicular to the plane containing the local field and the radial direction.
    e_perp2 = np.cross(e_parallel, ref_vec_radial)
    e_perp2 = e_perp2 / np.linalg.norm(e_perp2)

    # e_perp1 completes the right-handed system (e_perp1, e_perp2, e_parallel).
    e_perp1 = np.cross(e_perp2, e_parallel)

    # Step 5: Output the results, explaining their significance.
    # The prompt asks to "output each number in the final equation". This is
    # interpreted as printing the components of the resulting basis vectors.
    print("This script shows how to find the correct coordinate system for helicity analysis.")
    print("Your intuition is correct: the plane perpendicular to the *local mean magnetic field* should be used.")
    print("-" * 70)
    print(f"Assumed local mean magnetic field B_mean (nT):")
    print(f"Bx = {B_mean[0]:.4f}, By = {B_mean[1]:.4f}, Bz = {B_mean[2]:.4f}")
    print("\nThis field is not purely radial, as it has significant By and Bz components.")
    print("-" * 70)
    print("The basis vectors for the proper Field-Aligned Coordinate (FAC) system are:")

    print("\n1. The axis PARALLEL to the local magnetic field:")
    print(f"   e_parallel = [{e_parallel[0]:.4f}, {e_parallel[1]:.4f}, {e_parallel[2]:.4f}]")

    print("\n2. The plane PERPENDICULAR to the local magnetic field is defined by:")
    print(f"   e_perp1    = [{e_perp1[0]:.4f}, {e_perp1[1]:.4f}, {e_perp1[2]:.4f}]")
    print(f"   e_perp2    = [{e_perp2[0]:.4f}, {e_perp2[1]:.4f}, {e_perp2[2]:.4f}]")
    print("-" * 70)
    print("\nCONCLUSION:")
    print("To correctly calculate helicity, one must use the fluctuating magnetic field")
    print("components projected onto the 'e_perp1' and 'e_perp2' axes, not the original")
    print("By and Bz components (unless the field is perfectly radial).")

if __name__ == '__main__':
    calculate_field_aligned_coordinates()