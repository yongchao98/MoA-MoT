import textwrap

def identify_compound_a():
    """
    This function analyzes the provided reaction and NMR data to identify Compound A.
    """

    # --- Step 1: Define the provided data ---
    h_nmr = {
        8.69: "t, J = 5.7 Hz, 1H (likely -NH-CH2)",
        8.24: "s, 1H (likely aromatic CH or NH)",
        8.11: "s, 1H (likely aromatic CH or NH)",
        "7.37-7.22": "m, 5H (C6H5-, monosubstituted phenyl group)",
        4.73: "d, J = 6.0 Hz, 2H (likely -CH2-NH)",
        1.70: "s, 9H (tert-butyl group, -C(CH3)3)"
    }

    c_nmr_shifts = [156.89, 154.96, 152.80, 139.82, 130.16, 128.82, 127.85, 127.35, 102.23, 59.79, 43.52, 29.25]

    # --- Step 2: Analyze the reaction and NMR fragments ---
    analysis_text = """
    Step-by-Step Analysis:

    1. Reaction Analysis:
       - The reaction is a two-step nucleophilic substitution.
       - Step 1 uses tert-butyl hydrazine (H2N-NH-tBu) as a nucleophile.
       - Step 2 uses benzylamine (PhCH2-NH2) as a nucleophile.
       - This implies that the starting material, Compound A, has two leaving groups that are replaced sequentially.

    2. NMR Fragment Assignment:
       Based on the NMR data, we can identify the fragments incorporated into the final product.

       - Benzylamino group (from benzylamine):
         - 1H NMR: The 5H multiplet at 7.37-7.22 ppm is characteristic of a C6H5 group. The doublet at 4.73 ppm (2H) coupled to the triplet at 8.69 ppm (1H) confirms the -CH2-NH- moiety.
         - 13C NMR: Signals at 139.82, 128.82, 127.85, 127.35 (phenyl ring) and 43.52 (CH2) correspond to this group.

       - Tert-butyl group (from tert-butyl hydrazine):
         - 1H NMR: The 9H singlet at 1.70 ppm is the unmistakable signal of a tert-butyl group.
         - 13C NMR: Signals at 59.79 (quaternary C) and 29.25 (methyl C's) correspond to this group.

    3. Core Structure Determination:
       - Let's count the assigned 13C signals: 4 (phenyl) + 1 (CH2) + 1 (tBu quat C) + 1 (tBu methyl C) = 7 signals.
       - Total 13C signals given = 12.
       - Number of signals for the core skeleton = 12 - 7 = 5 signals.
       - The shifts of these 5 core carbons are: {core_shifts}.
       - The 1H NMR has two unassigned singlets (8.24 and 8.11 ppm), suggesting the core has two CH protons which are not adjacent. The other three core carbons must be quaternary.
       - A pyrimidine ring only has 4 carbons, which does not fit. A purine ring system has a 5-carbon framework (C2, C4, C5, C6, C8). This fits the data.

    4. Identification of Compound A:
       - A purine core with two leaving groups is the most plausible structure for Compound A.
       - A C2,C6-disubstituted purine product would have one CH group (at C8), which would be a singlet. The second singlet observed in the 1H NMR (8.24 or 8.11 ppm) can be attributed to an N-H proton on the purine ring (e.g., N9-H), which often appears as a sharp singlet in DMSO solvent.
       - The most common starting material for this type of synthesis is a dichloropurine. The reaction sequence is a classic double nucleophilic aromatic substitution on this scaffold.
    """

    core_carbon_shifts = sorted([s for s in c_nmr_shifts if s not in [139.82, 128.82, 127.85, 127.35, 59.79, 43.52, 29.25]], reverse=True)
    print(textwrap.dedent(analysis_text.format(core_shifts=core_carbon_shifts)).strip())
    
    # --- Step 3: State the conclusion ---
    conclusion = "Based on this analysis, the starting material, Compound A, is 2,6-Dichloropurine."
    print("\n" + "="*40)
    print("Conclusion")
    print("="*40)
    print(conclusion)

# Execute the analysis
identify_compound_a()