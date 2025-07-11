def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum by assigning its peaks to the functional groups
    of the most likely candidate molecule, C-L (Lidocaine analog).
    """
    
    # Structure C-L: N-(2-methylphenyl)-2-(diethylamino)acetamide
    
    print("Step-by-step analysis of the NMR spectrum corresponding to structure C-L:")
    print("=======================================================================")
    
    print("1. Peak at ~8.9 ppm (broad singlet, 1H):")
    print("   - This peak is assigned to the single proton of the amide group (-NH-).")
    print("   - Equation: Amide -N(H)-")
    print("     1H")
    print("-" * 30)

    print("2. Peak at ~7.1 ppm (multiplet, 4H):")
    print("   - These are the four protons on the substituted benzene ring.")
    print("   - Equation: Aromatic C-H")
    print("     4H")
    print("-" * 30)

    print("3. Peak at ~3.4 ppm (singlet, 2H):")
    print("   - This signal corresponds to the methylene (-CH2-) protons located between the")
    print("     carbonyl group (C=O) and the diethylamino group. They appear as a singlet")
    print("     because they have no adjacent protons.")
    print("   - Equation: -CO-CH2-N")
    print("                2H")
    print("-" * 30)
    
    print("4. Peak at ~2.8 ppm (quartet, 4H):")
    print("   - This quartet arises from the two methylene (-CH2-) groups of the two ethyl")
    print("     substituents. Each is split into a quartet by its 3 neighboring methyl protons.")
    print("   - Equation: -N-(CH2CH3)2")
    print("                  4H")
    print("-" * 30)
    
    print("5. Peak at ~2.3 ppm (singlet, 3H):")
    print("   - This singlet is from the methyl (-CH3) group attached directly to the aromatic ring.")
    print("   - Equation: Aromatic -CH3")
    print("                  3H")
    print("-" * 30)
    
    print("6. Peak at ~1.2 ppm (triplet, 6H):")
    print("   - This triplet corresponds to the six protons of the two methyl (-CH3) groups")
    print("     of the ethyl substituents. Each is split into a triplet by its 2 neighboring")
    print("     methylene protons.")
    print("   - Equation: -N-(CH2CH3)2")
    print("                        6H")
    print("-" * 30)

    print("\nConclusion: All peaks in the spectrum are accounted for by structure C-L.")

analyze_nmr_spectrum()