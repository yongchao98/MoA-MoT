import collections

def identify_compound_a():
    """
    This script identifies a starting material (Compound A) by analyzing the
    reaction pathway and the NMR spectral data of the final product.
    """

    # --- Step 1: Analyze data from the problem description ---
    
    # Product fragments based on reagents
    reagent_fragments = {
        'tert_butyl_hydrazine': '(CH3)3C-NH-NH-',
        'benzylamine': 'C6H5-CH2-NH-'
    }

    # ¹H NMR data (integrations)
    h_nmr_integrations = [1, 1, 1, 5, 2, 9]

    # ¹³C NMR data (number of signals)
    c_nmr_signals_total = 12

    print("Step 1: Analyzing product structure from reagents and NMR data.")
    print("-" * 60)

    # --- Step 2: Determine the number of C-H protons on the core structure ---

    # Protons from benzylamine fragment: 5 (aromatic) + 2 (CH2) + 1 (NH)
    protons_benzylamine = 5 + 2 + 1
    # Protons from tert-butyl hydrazine fragment: 9 (tert-butyl) + 2 (NH)
    protons_tBu_hydrazine = 9 + 1 + 1
    
    total_protons_from_fragments = protons_benzylamine + protons_tBu_hydrazine
    total_protons_observed = sum(h_nmr_integrations)

    core_protons = total_protons_observed - total_protons_from_fragments

    print(f"The ¹H NMR spectrum shows signals for a total of {total_protons_observed} protons ({' + '.join(map(str, h_nmr_integrations))}).")
    print(f"The known fragments from benzylamine and tert-butyl hydrazine account for {total_protons_from_fragments} protons.")
    print(f"Calculation: {total_protons_observed} (observed) - {total_protons_from_fragments} (fragments) = {core_protons} core protons.")
    print("Conclusion: The core of the molecule has no C-H bonds.\n")

    # --- Step 3: Determine the number of carbons in the core structure ---
    
    # Carbons from benzylamine fragment: 4 (unique aromatic C) + 1 (CH2)
    carbons_benzylamine = 4 + 1
    # Carbons from tert-butyl fragment: 1 (quaternary C) + 1 (methyl C)
    carbons_tBu = 1 + 1
    
    total_carbons_from_fragments = carbons_benzylamine + carbons_tBu
    core_carbons = c_nmr_signals_total - total_carbons_from_fragments
    
    print(f"The ¹³C NMR spectrum shows a total of {c_nmr_signals_total} unique carbon signals.")
    print(f"The known fragments account for {total_carbons_from_fragments} of these signals.")
    print(f"Calculation: {c_nmr_signals_total} (observed) - {total_carbons_from_fragments} (fragments) = {core_carbons} core carbons.")
    print("Conclusion: The core of the molecule is composed of 5 carbon atoms.\n")
    
    # --- Step 4: Deduce the identity of Compound A ---
    
    print("Step 2: Deducing the identity of the starting material (Compound A).")
    print("-" * 60)
    print("The starting material must satisfy the following criteria:")
    print(f"  - Has a {core_carbons}-carbon core.")
    print(f"  - Has {core_protons} C-H bonds on the core.")
    print("  - Has at least two reactive sites (leaving groups) for the sequential substitution.")
    
    print("\nA purine is a common heterocyclic structure with a 5-carbon, 4-nitrogen core.")
    print("While simple dichloropurines (like 2,6-dichloropurine) have a C-H bond, a trichlorinated version would not.")
    print("2,6,8-trichloropurine has a 5-carbon core, no C-H bonds, and three chloro leaving groups.")
    print("It is a common starting material for sequential nucleophilic substitutions at the 6, 2, and 8 positions.")
    print("\nTherefore, the most plausible identity for Compound A is:")
    
    final_answer = "2,6,8-trichloropurine"
    print(f"\nFinal Answer: {final_answer}")

if __name__ == '__main__':
    identify_compound_a()
