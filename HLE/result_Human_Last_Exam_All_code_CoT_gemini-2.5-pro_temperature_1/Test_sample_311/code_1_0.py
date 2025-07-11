def solve_synthesis_questions():
    """
    Solves questions about atom conservation in a chemical synthesis.
    """

    # Question 1: How many carbons from compound 11 are present in compound 1?
    # Step 1: Compound 11 is cis-2-butene-1,4-diol. It has a 4-carbon chain (C4).
    # Step 2: Ozonolysis of the double bond cleaves the 4-carbon chain into two identical 2-carbon fragments.
    # Step 3: The synthesis proceeds using one of these 2-carbon fragments (aldehyde 12).
    # Step 4: The subsequent reactions modify and elongate other parts of the molecule,
    # but the core 2-carbon fragment from compound 11 is not broken or lost. It is carried through to the final product, compound 1.
    carbons_from_11_in_1 = 2

    # Question 2: How many oxygens from compound 11 are present in compound 14?
    # Step 1: Compound 11 (HO-CH2-CH=CH-CH2-OH) has 2 oxygen atoms.
    # Step 2: After cleavage by ozonolysis, each 2-carbon fragment contains one of the original hydroxy groups (-OH).
    # Thus, the fragment used for the synthesis contains 1 oxygen atom from compound 11.
    # Step 3: This oxygen is protected (TESO-) and remains in the molecule through the synthesis of compounds 12, 10, 7, 13, and 14.
    oxygens_from_11_in_14 = 1

    # Question 3: How many nitrogens from compound 7 are present in compound 10?
    # Step 1: The question is chronologically reversed. Compound 10 is a precursor to 7. We determine the number of nitrogens conserved between them.
    # Step 2: Compound 10 is a nitroolefin, formed by reacting aldehyde 12 with MeNO2 (nitromethane). It contains one nitro group (-NO2) and therefore 1 nitrogen atom.
    # Step 3: Compound 7 is formed by an addition reaction to compound 10. The text states that the nitro group is later converted to an aldehyde via a Nef reaction (to make 14 from 13), which means the nitro group is present in compound 7.
    # Step 4: Therefore, both compound 10 and compound 7 contain 1 nitrogen atom, which originates from MeNO2.
    nitrogens_in_10_and_7 = 1
    
    # The final answer consists of the three numbers separated by commas.
    print(f"{carbons_from_11_in_1}, {oxygens_from_11_in_14}, {nitrogens_in_10_and_7}")

solve_synthesis_questions()