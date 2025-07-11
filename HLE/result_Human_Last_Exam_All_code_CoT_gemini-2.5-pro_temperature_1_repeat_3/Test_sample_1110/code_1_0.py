def analyze_boronic_acid_synthesis(eq_tmb):
    """
    Calculates the molar ratio of product to byproduct in a boronic acid synthesis.

    Args:
        eq_tmb (float): The equivalents of trimethyl borate (TMB) used.
    """
    # We can normalize the calculation by assuming we start with 1 mole of the aryl halide.
    moles_aryl_halide = 1.0

    # The reaction is 1:1, so 1 mole of aryl halide forms 1 mole of product.
    moles_product = moles_aryl_halide
    moles_tmb_reacted = moles_aryl_halide

    # Calculate total moles of TMB used and the excess amount.
    moles_tmb_total = moles_aryl_halide * eq_tmb
    moles_tmb_excess = moles_tmb_total - moles_tmb_reacted

    # During workup, excess TMB hydrolyzes to form boric acid byproduct B(OH)3.
    moles_boric_acid_byproduct = moles_tmb_excess

    print(f"Analysis for {eq_tmb} equivalents of trimethyl borate:")
    print(f"For every {moles_product} mole of desired boronic acid produced...")
    # The final equation shows the molar ratio of product to byproduct.
    print(f"Final Molar Ratio => Product : Byproduct")
    print(f"                      {moles_product:.1f}  :  {moles_boric_acid_byproduct:.1f}")
    print("-" * 30)

# --- Scenario from the problem description ---
print("Original reaction conditions resulted in two B-NMR signals.")
analyze_boronic_acid_synthesis(5.0)
print("The large amount of boric acid byproduct (4 moles per mole of product) is the likely cause.\n")


# --- Proposed solution (Choice D) ---
print("Proposed solution: Use a smaller excess of trimethyl borate.")
analyze_boronic_acid_synthesis(1.2)
print("This significantly reduces the byproduct, simplifying purification and likely solving the problem.")
