import math

def calculate_steric_virial_coefficient():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a standard monoclonal antibody (mAb) using a hard-sphere model.
    """
    # 1. Define the physical constants and standard parameters for a mAb.
    # Molar Mass (M) of a typical mAb in g/mol
    M = 150000.0
    # Hydrodynamic Radius (R) of a typical mAb in nm
    R_nm = 5.5
    # Convert Radius from nm to cm for unit consistency (1 nm = 1e-7 cm)
    R_cm = R_nm * 1e-7
    # Avogadro's number (N_A) in mol^-1
    N_A = 6.022e23
    # Pi
    pi = math.pi

    # 2. Calculate the steric-only B22 using the hard-sphere formula.
    # The formula for B22 in units of volume/mass (e.g., mL/g or cm^3/g) is:
    # B22_steric = (16 * pi * N_A * R^3) / (3 * M)
    numerator = 16 * pi * N_A * (R_cm**3)
    denominator = 3 * M
    B22_steric = numerator / denominator

    # 3. Print the parameters and the final result.
    # Since 1 cm^3 = 1 mL, the calculated value is in mL/g.
    print("Calculation of the Steric-Only Second Osmotic Virial Coefficient (B22, steric)")
    print("-" * 75)
    print("The calculation uses a hard-sphere model with standard parameters for a mAb.")
    
    print("\nParameters used in the final equation:")
    print(f"  Molar Mass (M)             = {M} g/mol")
    print(f"  Avogadro's Number (N_A)    = {N_A:.4g} /mol")
    print(f"  Hydrodynamic Radius (R)    = {R_nm} nm = {R_cm:.1e} cm")
    
    print("\nFinal Equation:")
    print(f"  B22, steric = (16 * pi * N_A * R^3) / (3 * M)")
    print(f"  B22, steric = (16 * {pi:.5f} * {N_A:.4g} * ({R_cm:.1e} cm)^3) / (3 * {M} g/mol)")

    print("\nResult:")
    print(f"  The second osmotic virial coefficient from steric-only behavior is: {B22_steric:.2f} mL/g")
    print("-" * 75)
    print("\nThis positive value represents the repulsive contribution from the protein's excluded")
    print("volume. The overall negative experimental value (-7.585 mL/g) indicates that")
    print("strong attractive interactions overcome this inherent repulsion at pH 5.")

calculate_steric_virial_coefficient()
<<<11.19>>>