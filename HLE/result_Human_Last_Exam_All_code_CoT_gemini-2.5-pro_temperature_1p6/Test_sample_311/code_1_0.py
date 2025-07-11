def solve_synthesis_questions():
    """
    This function solves the three questions based on the provided chemical synthesis description.
    1. How many carbons from compound 11 are present in compound 1?
    2. How many oxygens from compound 11 are present in compound 14?
    3. How many nitrogens from compound 7 are present in compound 10?
    """

    # Answer to question 1
    # Compound 11 (C4) is cleaved by ozonolysis into two C2 fragments (compound 12).
    # The synthesis uses one C2 fragment. No subsequent reactions break the C-C bond of this fragment.
    # So, 2 carbons from 11 are in the final product 1.
    carbons_11_in_1 = 2

    # Answer to question 2
    # Compound 11 (a diol) is cleaved into two fragments, each containing one of the original oxygen atoms.
    # The synthesis follows one fragment. This oxygen (in a protected form) is retained through the synthesis of 14.
    # So, 1 oxygen from 11 is in compound 14.
    oxygens_11_in_14 = 1

    # Answer to question 3
    # The question is chronologically backward (10 -> 7). We interpret it as "How many nitrogens in compound 10?".
    # Compound 10 is made from compound 12 and MeNO2. MeNO2 has one nitrogen atom.
    # So, compound 10 has 1 nitrogen.
    nitrogens_in_10 = 1

    # Print the answers separated by commas
    print(f"{carbons_11_in_1},{oxygens_11_in_14},{nitrogens_in_10}")

solve_synthesis_questions()
# <<<2,1,1>>>