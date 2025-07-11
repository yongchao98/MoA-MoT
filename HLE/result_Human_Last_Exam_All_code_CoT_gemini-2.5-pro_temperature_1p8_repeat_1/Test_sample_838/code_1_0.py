import textwrap

def explain_reaction_sequence():
    """
    This function explains the three-step chemical synthesis and identifies the final product.
    """
    
    # Wrap long lines for better readability
    wrapper = textwrap.TextWrapper(width=80)

    # Explanation of Step 1
    print("Step 1: E2 Elimination to form Product A")
    print("----------------------------------------")
    explanation1 = "The starting material, [(3S)-3-bromobutyl]benzene (IUPAC name: (S)-2-bromo-4-phenylbutane), is a secondary alkyl halide. It is treated with potassium tert-butoxide (t-BuOK), a strong, sterically bulky base, in a 60/40 cyclohexane/diethyl ether mixture. These conditions favor an E2 elimination reaction. Due to the large size of the base, it preferentially removes a proton from the more accessible, less sterically hindered position (the terminal methyl group), following the Hofmann rule."
    print("\n".join(wrapper.wrap(text=explanation1)))
    print("Product A: The major product is the less substituted alkene, 4-phenylbut-1-ene. The chiral center is destroyed in this step, so Product A is achiral.\n")

    # Explanation of Step 2
    print("Step 2: Hydroboration-Oxidation to form Product B")
    print("-------------------------------------------------")
    explanation2 = "Product A (4-phenylbut-1-ene) is subjected to hydroboration-oxidation. This two-step procedure involves an initial reaction with borane in THF (BH3·THF), followed by oxidation with hydrogen peroxide (H2O2) and sodium hydroxide (NaOH). This reaction results in the anti-Markovnikov addition of water across the double bond."
    print("\n".join(wrapper.wrap(text=explanation2)))
    print("Product B: The hydroxyl (-OH) group adds to the terminal, less-substituted carbon of the former double bond, yielding the primary alcohol, 4-phenylbutan-1-ol.\n")

    # Explanation of Step 3
    print("Step 3: Bromination to form Product C")
    print("--------------------------------------")
    explanation3 = "Product B (4-phenylbutan-1-ol) is treated with phosphorous tribromide (PBr3). This is a standard reagent that converts primary alcohols into primary alkyl bromides through a reaction that proceeds with an SN2 mechanism."
    print("\n".join(wrapper.wrap(text=explanation3)))
    print("Product C: The -OH group is replaced by a -Br atom to yield the final product, 1-bromo-4-phenylbutane.\n")

    # Final Product Identity and Chirality
    print("Final Product C: Identity and Chirality")
    print("----------------------------------------")
    print("IUPAC Name: 1-bromo-4-phenylbutane")
    
    chirality_explanation = "Chirality: The final product, 1-bromo-4-phenylbutane, is achiral. Although the starting material was chiral, the first step of the synthesis was an elimination reaction that formed a flat C=C double bond, which destroyed the molecule's only chiral center. Since no new chiral centers were generated in the subsequent steps, the final product is not chiral."
    print("\n".join(wrapper.wrap(text=chirality_explanation)))

explain_reaction_sequence()
