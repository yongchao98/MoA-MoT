def print_reaction_stoichiometry():
    """
    This function outlines the balanced chemical equation for the ethylation of
    2-Methyl-1,4-naphthalenediol and prints the stoichiometric coefficients.
    The reaction proceeds in two conceptual steps:
    1. Deprotonation by NaH: C11H10O2 + 2 NaH -> C11H8O2(Na)2 + 2 H2
    2. Alkylation by EtBr: C11H8O2(Na)2 + 2 C2H5Br -> C15H18O2 + 2 NaBr
    The overall balanced equation is derived from these steps.
    """
    
    # Define the stoichiometric coefficients for the overall balanced reaction:
    # C11H10O2 + 2 NaH + 2 C2H5Br -> C15H18O2 + 2 NaBr + 2 H2
    
    coeff_diol = 1
    coeff_NaH = 2
    coeff_EtBr = 2
    coeff_product = 1
    coeff_NaBr = 2
    coeff_H2 = 2
    
    print("The balanced chemical equation for the desired reaction is:")
    # The chemical formulas represent:
    # C11H10O2: 2-Methyl-1,4-naphthalenediol (Starting Material)
    # NaH: Sodium Hydride (Base)
    # C2H5Br: Ethyl Bromide (Alkylating Agent)
    # C15H18O2: 1,4-Diethoxy-2-methylnaphthalene (Product)
    # NaBr: Sodium Bromide (Byproduct)
    # H2: Hydrogen Gas (Byproduct)
    print(f"{coeff_diol} C11H10O2 + {coeff_NaH} NaH + {coeff_EtBr} C2H5Br -> {coeff_product} C15H18O2 + {coeff_NaBr} NaBr + {coeff_H2} H2")
    
    print("\nHere are the integer coefficients (the numbers) for the final balanced equation:")
    print(f"Coefficient for Starting Material (C11H10O2): {coeff_diol}")
    print(f"Coefficient for Sodium Hydride (NaH): {coeff_NaH}")
    print(f"Coefficient for Ethyl Bromide (C2H5Br): {coeff_EtBr}")
    print(f"Coefficient for Product (C15H18O2): {coeff_product}")
    print(f"Coefficient for Sodium Bromide (NaBr): {coeff_NaBr}")
    print(f"Coefficient for Hydrogen Gas (H2): {coeff_H2}")

print_reaction_stoichiometry()