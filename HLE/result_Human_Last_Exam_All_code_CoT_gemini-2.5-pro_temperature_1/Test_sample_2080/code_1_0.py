import math

def calculate_steric_virial_coefficient():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody (mAb) using a hard-sphere model.
    """
    # Plan:
    # 1. State the assumptions about mAb properties (M and Rh).
    # 2. Define constants and variables.
    # 3. Calculate B22,steric using the hard-sphere model formula.
    # 4. Print the final equation with all the numbers and the result.

    # Constants
    N_A = 6.022e23  # Avogadro's number, in mol^-1
    pi = math.pi

    # Assumed typical properties for a monoclonal antibody (mAb)
    M = 150000.0      # Molar Mass, in g/mol
    Rh_nm = 5.0       # Hydrodynamic Radius, in nm
    
    # Convert Hydrodynamic Radius from nm to cm for the calculation
    # (1 nm = 1e-7 cm)
    Rh_cm = Rh_nm * 1e-7

    # Calculate B22,steric using the formula for hard-sphere interactions.
    # The result will be in cm^3/g, which is equivalent to mL/g.
    numerator = 16 * pi * N_A * (Rh_cm**3)
    denominator = 3 * M
    B22_steric = numerator / denominator

    print("The second osmotic virial coefficient from Steric-only behavior (B22,steric) is calculated using a hard-sphere model.")
    print("Using typical properties for a monoclonal antibody: M = 150,000 g/mol and hydrodynamic radius (Rh) = 5.0 nm.")
    print("\nThe equation with the values plugged in is:")
    
    # Output the final equation with all the numbers and the result.
    # The format shows the formula, the values used, and the final answer.
    print(f"B22,steric = (16 * {pi:.4f} * {N_A:.3e} * ({Rh_cm:.1e} cm)^3) / (3 * {M} g/mol) = {B22_steric:.3f} mL/g")

calculate_steric_virial_coefficient()