import math

def calculate_steric_virial_coefficient():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody (mAb).
    """
    # Step 1: Define constants and typical mAb parameters.
    # Avogadro's number (mol^-1)
    N_A = 6.022e23
    # Typical molecular weight of an mAb (g/mol)
    M_w = 150000.0
    # Typical hard-sphere radius of an mAb (nm)
    R_hs_nm = 5.0
    # Convert radius to cm for volume calculation
    R_hs_cm = R_hs_nm * 1e-7

    # Step 2: Calculate the volume of a single mAb molecule (cm^3).
    V_molecule_cm3 = (4/3) * math.pi * (R_hs_cm**3)

    # Step 3: Calculate the specific volume (v_sp) of the hydrated mAb.
    # The result will be in cm^3/g, which is equivalent to mL/g.
    v_sp = (N_A * V_molecule_cm3) / M_w

    # Step 4: Calculate the steric-only second virial coefficient (B22_steric).
    # The formula is B22_steric = 4 * v_sp.
    B22_steric = 4 * v_sp

    # Step 5: Print the results, showing the numbers in the final equation.
    print(f"Assuming a typical mAb with M_w = {M_w:.0f} g/mol and R_hs = {R_hs_nm:.1f} nm.")
    print(f"The calculated specific volume (v_sp) is: {v_sp:.3f} mL/g")
    print("\nThe steric-only second osmotic virial coefficient (B22_steric) is calculated as:")
    print(f"B22_steric = 4 * v_sp")
    print(f"B22_steric = 4 * {v_sp:.3f} mL/g = {B22_steric:.3f} mL/g")

calculate_steric_virial_coefficient()