import numpy as np

def calculate_field_aligned_components():
    """
    Demonstrates the transformation of magnetic field data from a standard
    coordinate system (like RTN) to a Field-Aligned Coordinate (FAC) system.
    """
    # Let's assume a background magnetic field vector B at L1, given in a
    # standard coordinate system like RTN (Radial, Tangential, Normal).
    # A typical Parker spiral at 1 AU might give Bx ~ -4 nT, By ~ -4 nT.
    # B = [Br, Bt, Bn] in nT
    B_background_rtn = np.array([-4.0, -4.0, 1.5])

    # Let's assume we have a small magnetic field fluctuation (the AIC wave)
    # measured at the same time in the same RTN system.
    # b_wave = [br, bt, bn] in nT
    b_wave_rtn = np.array([0.1, -0.2, 0.3])

    print("--- Initial Data in RTN System ---")
    print(f"Background Magnetic Field B = {B_background_rtn} nT")
    print(f"Wave Fluctuation b = {b_wave_rtn} nT\n")

    # Step 1: Define the axes of the Field-Aligned Coordinate (FAC) system.
    # The new z-axis (e_parallel) is parallel to the background magnetic field.
    B_magnitude = np.linalg.norm(B_background_rtn)
    e_parallel = B_background_rtn / B_magnitude

    # The new y-axis (e_perp1) is perpendicular to both B and the RTN Normal axis ([0,0,1]).
    # This choice is a common convention to fix the orientation.
    e_perp1 = np.cross(e_parallel, np.array([0., 0., 1.]))
    e_perp1 = e_perp1 / np.linalg.norm(e_perp1)

    # The new x-axis (e_perp2) completes the right-handed system (e_perp2 = e_perp1 x e_parallel).
    e_perp2 = np.cross(e_perp1, e_parallel)

    print("--- Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    print(f"e_parallel (new z-axis) = {np.round(e_parallel, 4)}")
    print(f"e_perp1    (new y-axis) = {np.round(e_perp1, 4)}")
    print(f"e_perp2    (new x-axis) = {np.round(e_perp2, 4)}\n")


    # Step 2: Create the transformation matrix from RTN to FAC.
    # The rows of this matrix are the new basis vectors.
    # B_fac = T * B_rtn
    T = np.array([e_perp2, e_perp1, e_parallel])

    # Step 3: Transform the wave fluctuation vector into the FAC system.
    b_wave_fac = T.dot(b_wave_rtn)
    
    # Unpack the components for clarity
    b_perp2 = b_wave_fac[0]
    b_perp1 = b_wave_fac[1]
    b_parallel_comp = b_wave_fac[2]


    print("--- Transformed Wave Components in FAC System ---")
    print(f"Original wave vector b_rtn: [br={b_wave_rtn[0]}, bt={b_wave_rtn[1]}, bn={b_wave_rtn[2]}]")
    print(f"Transformed wave vector b_fac: [b_perp2={b_perp2:.4f}, b_perp1={b_perp1:.4f}, b_parallel={b_parallel_comp:.4f}] nT\n")

    # Step 4: Justify the final result for helicity calculation.
    # The normalized magnetic helicity calculation requires the complex signal
    # B_perp(t) = b_perp1(t) + i*b_perp2(t).
    # The code below shows the numbers that would be used at one time step.
    print("--- Final Components for Helicity Calculation ---")
    print("To calculate the normalized magnetic helicity, you must use the components of the wave perpendicular to the background field.")
    print("For this specific wave fluctuation, the required components are:")
    print(f"Component 1 (b_perp1): {b_perp1:.4f}")
    print(f"Component 2 (b_perp2): {b_perp2:.4f}")
    print("\nThese are the values that replace the generic 'Y' and 'Z' components in the helicity formula for a rigorous calculation.")

# Execute the function
calculate_field_aligned_components()