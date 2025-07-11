import math

def calculate_steric_virial_coefficient():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody (mAb).
    """
    # Standard physical constants and typical mAb parameters
    # Avogadro's number (mol^-1)
    Na = 6.022e23
    # Typical molar mass of a mAb (g/mol)
    M = 150000.0
    # Typical hydrodynamic radius of a mAb (nm)
    Rh_nm = 5.0
    # Convert radius from nanometers to centimeters (1 nm = 1e-7 cm)
    Rh_cm = Rh_nm * 1e-7

    # The formula for the steric-only second osmotic virial coefficient (B22) is:
    # B22_steric = (16 * pi * Na * Rh^3) / (3 * M)
    # The result will be in cm³/g, which is equivalent to mL/g.
    b22_steric = (16 * math.pi * Na * (Rh_cm**3)) / (3 * M)

    # Output the equation with the specific numbers used
    print("The equation to calculate the steric-only second osmotic virial coefficient (B22_steric) with the values plugged in is:")
    # The f-string below prints the full expression used in the calculation.
    print(f"B22_steric = (16 * {math.pi} * {Na} * ({Rh_cm})**3) / (3 * {M})")
    
    # Print the final result
    # The problem context implies strong attractions, leading to a negative total B22.
    # The steric component, however, is always positive (repulsive).
    print(f"\nThe calculated second osmotic virial coefficient from steric-only behavior is {b22_steric:.3f} mL/g.")

# Run the calculation and print the results
calculate_steric_virial_coefficient()