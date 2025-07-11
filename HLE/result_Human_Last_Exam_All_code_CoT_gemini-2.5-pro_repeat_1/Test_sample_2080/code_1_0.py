import math

def calculate_b22_steric():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody using the hard-sphere model.
    """

    # --- Parameters ---
    # Molar Mass (M) of a typical mAb in g/mol
    M_g_per_mol = 150000.0

    # Hydrodynamic Radius (R_H) of a typical mAb in nm
    R_H_nm = 5.5
    # Convert R_H to cm for calculation (1 nm = 1e-7 cm)
    R_H_cm = R_H_nm * 1e-7

    # Avogadro's number (N_A) in mol^-1
    N_A = 6.022e23

    # Pi constant
    pi = math.pi

    # --- Calculation ---
    # The formula for the steric-only second virial coefficient (B22) in mL/g is:
    # B22_steric = (16 * pi * N_A * R_H^3) / (3 * M)
    # The result will be in cm^3/g, which is equivalent to mL/g.

    numerator = 16 * pi * N_A * (R_H_cm ** 3)
    denominator = 3 * M_g_per_mol
    b22_steric = numerator / denominator

    # --- Output ---
    print("This script calculates the second osmotic virial coefficient from steric-only behavior (B22_steric).")
    print("It uses a hard-sphere model with typical parameters for a monoclonal antibody.")
    print("\nParameters used:")
    print(f"Molar Mass (M): {M_g_per_mol} g/mol")
    print(f"Hydrodynamic Radius (R_H): {R_H_nm} nm ({R_H_cm:.1e} cm)")
    print(f"Avogadro's Number (N_A): {N_A:.3e} mol^-1")

    print("\nCalculation based on the formula: B22_steric = (16 * pi * N_A * R_H^3) / (3 * M)")
    print("Substituting the values into the equation:")
    # Using f-string formatting to display the numbers in the equation
    print(f"B22_steric = (16 * {pi:.5f} * {N_A:.3e} mol^-1 * ({R_H_cm:.1e} cm)^3) / (3 * {M_g_per_mol} g/mol)")

    print(f"\nCalculated B22_steric value: {b22_steric:.3f} mL/g")

if __name__ == "__main__":
    calculate_b22_steric()
