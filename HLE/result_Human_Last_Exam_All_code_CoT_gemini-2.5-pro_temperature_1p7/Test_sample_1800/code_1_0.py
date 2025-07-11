import math

def calculate_catalyst_composition():
    """
    Calculates the weight percentage of Ni in a Ni-Ceria catalyst
    for a given optimal atomic ratio based on scientific literature.
    """
    # --- Step 1: Define Constants ---
    # Atomic masses (g/mol) are fundamental constants for this calculation.
    ATOMIC_MASS_NI = 58.69
    ATOMIC_MASS_CE = 140.12

    # From literature, a highly effective atomic ratio (Ni/Ce) is around 0.15.
    # This value represents a balance between active sites and catalyst stability.
    optimal_atomic_ratio_ni_ce = 0.15

    print("--- Analysis of Optimal Ni/Ce Ratio in Ni-Ceria Catalysts ---")
    print(f"\nThe target optimal atomic ratio (Ni atoms / Ce atoms) is: {optimal_atomic_ratio_ni_ce}\n")

    # --- Step 2: Calculate Mass Components based on the Atomic Ratio ---
    # To convert an atomic ratio to a weight ratio, we imagine having a certain number of atoms.
    # Let's assume we have 1 mole of Cerium (Ce) atoms.
    moles_ce = 1.0
    # Based on the ratio, the moles of Nickel (Ni) would be:
    moles_ni = moles_ce * optimal_atomic_ratio_ni_ce

    # Now, calculate the mass of each component for these mole amounts.
    # Mass = moles * atomic mass
    mass_ni = moles_ni * ATOMIC_MASS_NI
    mass_ce = moles_ce * ATOMIC_MASS_CE
    
    print("--- Equation for Mass Calculation ---")
    print(f"Based on a composition of {moles_ni:.2f} moles of Ni and {moles_ce:.1f} mole of Ce:")
    print(f"Mass of Ni = {moles_ni:.2f} * {ATOMIC_MASS_NI} = {mass_ni:.2f} g")
    print(f"Mass of Ce = {moles_ce:.1f} * {ATOMIC_MASS_CE} = {mass_ce:.2f} g\n")

    # --- Step 3: Calculate Total Mass and Final Weight Percentage ---
    # The weight percentage (wt%) of Ni is a common way to describe catalyst composition in practice.
    # wt% Ni = (mass of Ni / (mass of Ni + mass of Ce)) * 100
    total_mass = mass_ni + mass_ce
    weight_percent_ni = (mass_ni / total_mass) * 100

    print("--- Final Equation for Ni Weight Percentage (wt%) ---")
    print(f"Total Mass = {mass_ni:.2f} (Mass Ni) + {mass_ce:.2f} (Mass Ce) = {total_mass:.2f} g")
    print(f"Ni wt% = ({mass_ni:.2f} / {total_mass:.2f}) * 100 = {weight_percent_ni:.2f} %\n")
    
    print(f"Conclusion: An optimal atomic Ni/Ce ratio of {optimal_atomic_ratio_ni_ce} corresponds to roughly {weight_percent_ni:.2f} wt% Nickel in the catalyst.")

# Execute the function
calculate_catalyst_composition()