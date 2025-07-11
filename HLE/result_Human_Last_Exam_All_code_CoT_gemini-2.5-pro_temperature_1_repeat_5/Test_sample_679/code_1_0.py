def solve_nmr_puzzle():
    """
    This function explains the step-by-step reasoning to identify the compound
    and prints the final IUPAC name.
    """
    
    formula = "C7H14"
    nmr_data = "145(s), 112(t), 48(t), 27(d), 22(q), 21(q)"

    print("Step 1: Analyze the molecular formula.")
    print(f"The formula is {formula}. This corresponds to one degree of unsaturation (a double bond or a ring).")
    print("-" * 20)

    print("Step 2: Analyze the 13C NMR data.")
    print(f"The signals at 145(s) and 112(t) indicate a C=C double bond.")
    print("The 145(s) signal is a quaternary carbon (C), and the 112(t) signal is a CH2 group.")
    print("This confirms the presence of a terminal alkene group: >C=CH2.")
    print("-" * 20)

    print("Step 3: Account for all atoms and identify fragments.")
    print("The formula has 7 carbons, but there are only 6 signals. This implies two carbons are chemically equivalent or overlapping.")
    print("Based on the multiplicities, the fragments are:")
    print("- One quaternary C (from 145s)")
    print("- One alkene CH2 (from 112t)")
    print("- One alkane CH2 (from 48t)")
    print("- One alkane CH (from 27d)")
    print("- Three CH3 groups (from the two q signals, with one signal representing two equivalent methyls)")
    print("These fragments correctly sum to C7H14.")
    print("-" * 20)

    print("Step 4: Assemble the structure.")
    print("The fragments can be assembled into an ethyl group (-CH2CH3) and an isopropyl group (-CH(CH3)2).")
    print("These two groups must be attached to the quaternary carbon of the >C=CH2 fragment.")
    print("Structure: CH2=C(CH2CH3)(CH(CH3)2)")
    print("-" * 20)
    
    print("Step 5: Determine the IUPAC name.")
    print("The longest carbon chain containing the double bond is a 4-carbon chain (but-1-ene).")
    print("To maximize the number of substituents, we choose the chain that gives an ethyl group at position 2 and a methyl group at position 3.")
    print("Numbering starts from the CH2= end to give the double bond the lowest number (1).")
    print("The substituents are 'ethyl' at C2 and 'methyl' at C3.")
    print("Alphabetizing the substituents gives '2-ethyl-3-methyl'.")
    final_name = "2-ethyl-3-methylbut-1-ene"
    print(f"The final IUPAC name is: {final_name}")
    
solve_nmr_puzzle()
<<<2-ethyl-3-methylbut-1-ene>>>