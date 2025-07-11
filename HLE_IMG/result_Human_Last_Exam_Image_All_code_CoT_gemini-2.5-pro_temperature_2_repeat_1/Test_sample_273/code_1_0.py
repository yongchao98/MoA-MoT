import textwrap

def determine_absolute_configuration():
    """
    This function explains the step-by-step process of determining the absolute
    configuration of the molecule shown in the image.
    """

    print("Step 1: Identify the Main Chain and Chiral Centers")
    explanation1 = textwrap.dedent("""\
        The molecule contains both an alcohol (-OH) and an amine (-NH2) group. According to IUPAC nomenclature rules, the alcohol has higher priority and defines the main chain and numbering. The longest carbon chain containing the C-OH group has five carbons. Numbering is started from the end closest to the -OH group.
        
        The correct IUPAC name is: 5-amino-3-ethyl-4-methylpentan-2-ol
        
        Based on this, we can identify three chiral centers (carbons with four different substituents) at positions 2, 3, and 4 of the pentane chain.
        C1H3 - C2(OH) - C3(Et) - C4(Me) - C5H2NH2
        """)
    print(explanation1)
    
    print("\nStep 2: Assign Priorities and Determine Configuration for Each Center")
    print("We will use the Cahn-Ingold-Prelog (CIP) rules to assign priorities (1=highest, 4=lowest) to the substituents on each chiral center. The 3D structure is interpreted from the 2D drawing.")
    
    print("\n--- Analysis of C2 ---")
    explanation_c2 = textwrap.dedent("""\
        The drawing at C2 is ambiguous, showing both an -OH and a -CH3 group with dashed bonds. This is a likely drawing error. We will proceed by interpreting the stereochemistry based on the standard convention for a backbone. The diagram shows the -OH group on a dash. This implies the implied H atom (the fourth substituent) must be on a wedge.
        - Priorities (1=highest): 1: -OH, 2: -C3(...), 3: -C1(CH3), 4: -H
        - Orientation: The lowest priority group (-H) is on a wedge (pointing towards the viewer).
        - Path 1 -> 2 -> 3: This traces a Counter-Clockwise (S) path.
        - Result: Since the lowest priority group points towards the viewer, we reverse the designation.
        """)
    print(explanation_c2)
    print("Configuration at C2 is (R).")

    print("\n--- Analysis of C3 ---")
    explanation_c3 = textwrap.dedent("""\
        The drawing shows the Ethyl group (-CH2CH3) on a wedge. This implies the H atom is on a dash.
        - Priorities (1=highest): 1: -C2(...), 2: -C4(...), 3: -Ethyl, 4: -H
        - Orientation: The lowest priority group (-H) is on a dash (pointing away from the viewer).
        - Path 1 -> 2 -> 3: This traces a Clockwise (R) path.
        - Result: Since the lowest priority group points away from the viewer, we keep the designation.
        """)
    print(explanation_c3)
    print("Configuration at C3 is (R).")

    print("\n--- Analysis of C4 ---")
    explanation_c4 = textwrap.dedent("""\
        The drawing shows the Methyl group (-CH3) on a dash. This implies the H atom is on a wedge.
        - Priorities (1=highest): 1: -CH2NH2, 2: -C3(...), 3: -CH3, 4: -H
        - Orientation: The lowest priority group (-H) is on a wedge (pointing towards the viewer).
        - Path 1 -> 2 -> 3: This traces a Counter-Clockwise (S) path.
        - Result: Since the lowest priority group points towards the viewer, we reverse the designation.
        """)
    print(explanation_c4)
    print("Configuration at C4 is (R).")

    print("\nStep 3: Final Absolute Configuration")
    print("Combining the results for each chiral center, the absolute configuration is:")
    print("The configuration at carbon 2 is R.")
    print("The configuration at carbon 3 is R.")
    print("The configuration at carbon 4 is R.")
    print("\nFinal Answer: (2R, 3R, 4R)")

determine_absolute_configuration()