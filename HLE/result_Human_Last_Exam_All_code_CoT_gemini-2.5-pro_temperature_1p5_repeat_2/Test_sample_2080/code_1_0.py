import numpy as np

def calculate_steric_b22():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody (mAb).
    """
    # --- Define constants and parameters ---
    # Avogadro's number (mol^-1)
    N_A = 6.022e23
    
    # Assumed molecular weight for a typical mAb in g/mol (150 kDa)
    M_w_g_mol = 150000.0
    
    # Assumed hydrodynamic radius for a typical mAb in nm
    R_h_nm = 5.5
    # Convert hydrodynamic radius to cm for consistent units (1 nm = 1e-7 cm)
    R_h_cm = R_h_nm * 1e-7

    # --- Calculate B22 steric ---
    # The formula for the steric contribution to B22 in mL/g is:
    # B22_steric = (16 * pi * N_A * R_h^3) / (3 * M_w)
    
    numerator = 16 * np.pi * N_A * (R_h_cm**3)
    denominator = 3 * M_w_g_mol
    
    b22_steric = numerator / denominator

    # --- Print the output ---
    print("This script calculates the second osmotic virial coefficient from steric-only behavior (B22_steric).")
    print("\nThe formula used is: B22_steric = (16 * pi * N_A * R_h^3) / (3 * M_w)")
    print("\nUsing typical values for a monoclonal antibody:")
    print(f"  - Avogadro's number (N_A): {N_A:.3e} mol^-1")
    print(f"  - Molecular Weight (M_w): {M_w_g_mol} g/mol")
    print(f"  - Hydrodynamic Radius (R_h): {R_h_nm} nm = {R_h_cm:.1e} cm")

    print("\nPlugging the values into the equation:")
    # Note: np.pi is used for pi in the calculation, but we print its approximate value for clarity.
    print(f"B22_steric = (16 * {np.pi:.4f} * {N_A:.3e} mol^-1 * ({R_h_cm:.1e} cm)^3) / (3 * {M_w_g_mol} g/mol)")

    print(f"\nResult:")
    print(f"The second osmotic virial coefficient from Steric-only behavior is {b22_steric:.3f} mL/g.")

if __name__ == "__main__":
    calculate_steric_b22()