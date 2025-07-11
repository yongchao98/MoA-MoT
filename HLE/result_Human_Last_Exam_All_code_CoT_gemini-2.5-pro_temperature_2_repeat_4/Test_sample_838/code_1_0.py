import textwrap

def solve_synthesis():
    """
    This script explains the step-by-step chemical synthesis and identifies the final product C.
    """

    # Explanation of Step 1
    print("--- Analysis of the 3-Step Synthesis ---")
    print("\nStep 1: [(3S)-3-bromobutyl]benzene to Product A")
    step1_explanation = """
    The starting material is [(3S)-3-bromobutyl]benzene. It is a secondary alkyl halide. The reagent is potassium tert-butoxide (t-BuOK), which is a sterically bulky, strong base. These conditions strongly favor an E2 elimination reaction over substitution.

    In an E2 elimination, the base removes a proton from a carbon adjacent (beta) to the carbon bearing the leaving group (bromine). There are two types of beta-protons: those on C2 and those on C4. Because t-BuOK is a bulky base, it preferentially removes the less sterically hindered proton. The protons on C4 (the terminal methyl group) are more accessible than the protons on C2.

    This leads to the formation of the 'Hofmann product' (the less substituted alkene). The double bond forms between C3 and C4. The product, A, is 4-phenylbut-1-ene. The original chiral center at C3 is eliminated in this reaction, so Product A is achiral.
    """
    print(textwrap.dedent(step1_explanation).strip())

    # Explanation of Step 2
    print("\nStep 2: Product A to Product B")
    step2_explanation = """
    Product A, 4-phenylbut-1-ene, is treated with borane (BH3) in THF, followed by oxidation with hydrogen peroxide (H2O2) and sodium hydroxide (NaOH). This is a hydroboration-oxidation reaction.

    This reaction adds water across the double bond with anti-Markovnikov regioselectivity. This means the hydroxyl (-OH) group adds to the less substituted carbon atom of the double bond (C1 of the but-1-ene). The hydrogen adds to the more substituted carbon (C2).

    The result is the formation of a primary alcohol. Product B is 4-phenylbutan-1-ol. This molecule is achiral.
    """
    print(textwrap.dedent(step2_explanation).strip())

    # Explanation of Step 3
    print("\nStep 3: Product B to Product C")
    step3_explanation = """
    Product B, 4-phenylbutan-1-ol, is treated with phosphorous tribromide (PBr3). This is a standard reagent used to convert primary and secondary alcohols into their corresponding alkyl bromides.

    The reaction proceeds via an SN2-type mechanism, where the hydroxyl group is converted into a good leaving group and then displaced by a bromide ion. The -OH group on C1 is replaced by a -Br atom.

    The final product, C, is 1-bromo-4-phenylbutane.
    """
    print(textwrap.dedent(step3_explanation).strip())

    # Final Answer
    print("\n--- Final Product Identification ---")
    final_answer = """
    The final product, C, is a primary alkyl bromide.

    IUPAC Name: 1-bromo-4-phenylbutane

    The numbers in the IUPAC name are 1 and 4.

    Chirality Explanation: The final product, 1-bromo-4-phenylbutane, is achiral. A molecule is chiral if it is non-superimposable on its mirror image, which typically arises from the presence of a chiral center (a carbon atom bonded to four different groups). In 1-bromo-4-phenylbutane, no such chiral center exists:
    - C1 is bonded to Br, H, H, and the rest of the alkyl chain.
    - C2 is bonded to C1, C3, H, and H.
    - C3 is bonded to C2, C4, H, and H.
    - C4 is bonded to C3, a phenyl group, H, and H.
    Since no carbon atom is bonded to four unique groups, the molecule is achiral.
    """
    print(textwrap.dedent(final_answer).strip())

if __name__ == "__main__":
    solve_synthesis()