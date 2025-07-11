import math

def calculate_steric_b22():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody.
    """
    # --- Parameters ---
    # Molar mass of a typical mAb in g/mol
    Mw = 150000.0
    # Hydrodynamic radius of a typical mAb in nm
    Rh_nm = 5.0
    # Avogadro's number in mol^-1
    NA = 6.022e23

    # --- Conversions ---
    # Convert radius from nanometers (nm) to centimeters (cm) for volume calculation
    # 1 nm = 1e-7 cm
    Rh_cm = Rh_nm * 1e-7

    # --- Calculation ---
    # Formula for the steric (hard-sphere) second osmotic virial coefficient (B22)
    # B22_steric = (16 * pi * Rh^3 * NA) / (3 * Mw)
    # The result will be in cm^3/g, which is equivalent to mL/g.
    
    numerator = 16 * math.pi * (Rh_cm**3) * NA
    denominator = 3 * Mw
    b22_steric = numerator / denominator

    # --- Output ---
    print("Calculating the second osmotic virial coefficient from Steric-only behavior (B22_steric):")
    print("-" * 80)
    print("Formula: B22_steric = (16 * pi * Rh^3 * NA) / (3 * Mw)")
    print("\nUsing typical values for a monoclonal antibody:")
    print(f"  Molar Mass (Mw) = {Mw:.1f} g/mol")
    print(f"  Hydrodynamic Radius (Rh) = {Rh_nm:.1f} nm = {Rh_cm:.1e} cm")
    print(f"  Avogadro's Number (NA) = {NA:.3e} mol^-1")
    print("-" * 80)
    
    # Print the final equation with all numbers
    print("Final Equation:")
    print(f"B22_steric = (16 * {math.pi:.4f} * ({Rh_cm:.1e} cm)^3 * {NA:.3e} mol^-1) / (3 * {Mw:.1f} g/mol)")
    
    # Print the final result
    print(f"\nCalculated B22_steric = {b22_steric:.3f} mL/g")

calculate_steric_b22()
<<<8.408>>>