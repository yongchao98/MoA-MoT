import sys

def calculate_reaction_stoichiometry():
    """
    Calculates reagent amounts and theoretical yield for the ethylation of
    2-Methyl-1,4-naphthalenediol.
    """

    # --- Molecular Weights (g/mol) ---
    mw_sm = 174.20      # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 23.998     # Reagent 1: Sodium Hydride (NaH)
    mw_etbr = 108.97    # Reagent 2: Ethyl Bromide (C2H5Br)
    mw_product = 230.31 # Product: 2-Methyl-1,4-diethoxynaphthalene (C15H18O2)

    # --- Reaction Parameters from User ---
    mass_sm = 10.0      # Mass of starting material in grams
    eq_nah_used = 2.5   # Equivalents of NaH the student used
    eq_etbr_used = 3.0  # Equivalents of EtBr the student used

    # --- Stoichiometric Calculations ---
    moles_sm = mass_sm / mw_sm

    # Stoichiometry from the balanced equation
    # 1 mol of SM requires 2 mol of base and 2 mol of ethylating agent.
    stoich_coeff_sm = 1
    stoich_coeff_nah = 2
    stoich_coeff_etbr = 2
    stoich_coeff_product = 1
    stoich_coeff_nabr = 2
    stoich_coeff_h2 = 2


    # --- Print Results ---
    print("--- Reaction Stoichiometry Analysis ---")
    print(f"Starting Material: 2-Methyl-1,4-naphthalenediol (MW: {mw_sm} g/mol)")
    print(f"Starting Mass: {mass_sm:.1f} g")
    print(f"Calculated Moles of Starting Material: {moles_sm:.4f} mol\n")

    print("Balanced Equation:")
    # The prompt requests to output each number in the final equation.
    print(f"{stoich_coeff_sm} C11H10O2 + {stoich_coeff_nah} NaH + {stoich_coeff_etbr} C2H5Br  ->  {stoich_coeff_product} C15H18O2 + {stoich_coeff_nabr} NaBr + {stoich_coeff_h2} H2\n")

    print("--- Reagent Calculation ---")
    # Sodium Hydride
    moles_nah_stoich = moles_sm * stoich_coeff_nah
    mass_nah_stoich = moles_nah_stoich * mw_nah
    moles_nah_planned = moles_sm * eq_nah_used
    mass_nah_planned = moles_nah_planned * mw_nah
    print(f"Sodium Hydride (NaH), MW: {mw_nah:.2f} g/mol")
    print(f"- Stoichiometric amount ({stoich_coeff_nah}.0 eq): {moles_nah_stoich:.4f} mol, or {mass_nah_stoich:.2f} g")
    print(f"- Planned amount ({eq_nah_used} eq): {moles_nah_planned:.4f} mol, or {mass_nah_planned:.2f} g\n")

    # Ethyl Bromide
    moles_etbr_stoich = moles_sm * stoich_coeff_etbr
    mass_etbr_stoich = moles_etbr_stoich * mw_etbr
    moles_etbr_planned = moles_sm * eq_etbr_used
    mass_etbr_planned = moles_etbr_planned * mw_etbr
    print(f"Ethyl Bromide (EtBr), MW: {mw_etbr:.2f} g/mol")
    print(f"- Stoichiometric amount ({stoich_coeff_etbr}.0 eq): {moles_etbr_stoich:.4f} mol, or {mass_etbr_stoich:.2f} g")
    print(f"- Planned amount ({eq_etbr_used} eq): {moles_etbr_planned:.4f} mol, or {mass_etbr_planned:.2f} g\n")


    print("--- Theoretical Yield ---")
    # Theoretical yield is based on the limiting reagent, which is the starting material.
    moles_product_theoretical = moles_sm * stoich_coeff_product
    mass_product_theoretical = moles_product_theoretical * mw_product
    print(f"Product: 2-Methyl-1,4-diethoxynaphthalene (MW: {mw_product:.2f} g/mol)")
    print(f"- Max theoretical yield: {moles_product_theoretical:.4f} mol, or {mass_product_theoretical:.2f} g")

if __name__ == '__main__':
    calculate_reaction_stoichiometry()