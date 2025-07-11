def analyze_butyllithium_stoichiometry():
    """
    Analyzes the effect of n-BuLi concentration errors on reaction stoichiometry.

    This script demonstrates why using a precise amount of n-BuLi is critical.
    It simulates a common scenario in a chemistry lab where the assumed molarity
    of a reagent differs from its actual molarity, leading to unintended side reactions.
    """

    # --- Properties of the Starting Material (SM) ---
    # SM is 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    mw_sm = 317.34  # g/mol

    # --- Experimental Parameters ---
    mass_sm_g = 5.0  # grams of starting material
    target_eq_nBuLi = 1.05  # The intended stoichiometry
    
    # Molarity of n-BuLi solution AS ASSUMED from the bottle's label
    assumed_nBuLi_M = 2.5  # mol/L
    
    # Molarity of n-BuLi solution AS IT ACTUALLY IS (e.g., a fresh bottle)
    # This value is unknown to the experimentalist but is the source of the error.
    actual_nBuLi_M = 2.7  # mol/L

    print("--- Reaction Plan Analysis ---")
    print(f"Starting Material (SM): 2-bromo-4-chloro-1-iodobenzene (MW = {mw_sm:.2f} g/mol)")
    print(f"Planned amount of SM: {mass_sm_g} g\n")

    # --- Calculations based on ASSUMED Molarity ---
    moles_sm = mass_sm_g / mw_sm
    target_moles_nBuLi = moles_sm * target_eq_nBuLi
    vol_to_add_mL = (target_moles_nBuLi / assumed_nBuLi_M) * 1000

    print("--- Calculation Based on Assumed Concentration ({:.2f} M) ---".format(assumed_nBuLi_M))
    print(f"1. Moles of starting material = {mass_sm_g} g / {mw_sm:.2f} g/mol = {moles_sm:.4f} mol")
    print(f"2. Target moles of n-BuLi ({target_eq_nBuLi} eq) = {moles_sm:.4f} mol * {target_eq_nBuLi} = {target_moles_nBuLi:.4f} mol")
    print(f"3. Volume of {assumed_nBuLi_M:.2f} M n-BuLi to add = {target_moles_nBuLi:.4f} mol / {assumed_nBuLi_M:.2f} M = {vol_to_add_mL/1000:.4f} L = {vol_to_add_mL:.2f} mL\n")

    # --- Revealing the ACTUAL amount added ---
    actual_moles_added = (vol_to_add_mL / 1000) * actual_nBuLi_M
    actual_eq_added = actual_moles_added / moles_sm

    print("--- Reality with Actual Concentration ({:.2f} M) ---".format(actual_nBuLi_M))
    print(f"The chemist adds {vol_to_add_mL:.2f} mL, believing they are adding {target_moles_nBuLi:.4f} moles.")
    print(f"1. Actual moles of n-BuLi added = {vol_to_add_mL/1000:.4f} L * {actual_nBuLi_M:.2f} M = {actual_moles_added:.4f} mol")
    print(f"2. Actual equivalents added = {actual_moles_added:.4f} mol / {moles_sm:.4f} mol = {actual_eq_added:.2f} eq\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Instead of the intended {target_eq_nBuLi} equivalents, {actual_eq_added:.2f} equivalents were actually added.")
    print("This significant excess of n-BuLi is highly reactive and causes side reactions,")
    print("leading to the formation of byproducts and the observed second boron signal in the NMR.")
    print("\nSOLUTION: Titrate the n-BuLi solution immediately before use to know its precise concentration.")
    print("This allows for the addition of the correct stoichiometric amount, preventing side reactions.")

if __name__ == '__main__':
    analyze_butyllithium_stoichiometry()

<<<C>>>