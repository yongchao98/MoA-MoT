import math

def calculate_b22_steric():
    """
    Calculates the steric-only second osmotic virial coefficient (B_22_steric)
    for a typical monoclonal antibody using the hard-sphere model.
    """

    # --- Constants and Parameters ---
    # Molar mass (M_w) of a typical mAb in g/mol
    Mw = 150000.0
    # Hydrodynamic radius (r_h) of a typical mAb in cm (5.5 nm = 5.5e-7 cm)
    r_h_cm = 5.5e-7
    # Avogadro's constant (N_A) in mol^-1
    NA = 6.022e23
    # Pi
    pi = math.pi

    # --- Calculation ---
    # The formula for the steric contribution to the second virial coefficient is:
    # B_22_steric = (16/3 * pi * r_h^3 * N_A) / M_w
    b22_steric = (16/3 * pi * (r_h_cm**3) * NA) / Mw

    # --- Output ---
    print("This script calculates the second osmotic virial coefficient from steric-only behavior (B_22^steric).")
    print("\nModel and assumptions:")
    print("- The monoclonal antibody (mAb) is treated as a hard sphere.")
    print("- Typical physical properties for an mAb are used.")
    
    print("\nParameters used:")
    print(f"- Molar Mass (M_w): {Mw} g/mol")
    print(f"- Hydrodynamic Radius (r_h): {r_h_cm * 1e7:.1f} nm = {r_h_cm} cm")
    print(f"- Avogadro's Constant (N_A): {NA:.4g} mol^-1")
    print(f"- Pi (π): {pi:.6f}")

    print("\nCalculation:")
    print(f"B_22^steric = (16/3 * {pi:.6f} * ({r_h_cm})^3 cm^3 * {NA:.4g} mol^-1) / {Mw} g/mol")
    
    numerator_val = 16/3 * pi * (r_h_cm**3) * NA
    print(f"B_22^steric = ({numerator_val:.4g} cm^3/mol) / {Mw} g/mol")
    
    print(f"\nThe calculated second osmotic virial coefficient from steric-only behavior is:")
    # The units cm^3/g are equivalent to mL/g.
    print(f"B_22^steric = {b22_steric:.2f} mL/g")
    
    # Context using the provided experimental data
    b22_experimental = -7.585
    b22_electrostatic = b22_experimental - b22_steric
    print("\nFor context:")
    print(f"The experimental B_22 (total) is {b22_experimental} mL/g.")
    print(f"This implies an electrostatic contribution (B_22_electrostatic = B_22_total - B_22_steric) of approximately {b22_electrostatic:.2f} mL/g, indicating strong net attraction.")

calculate_b22_steric()
