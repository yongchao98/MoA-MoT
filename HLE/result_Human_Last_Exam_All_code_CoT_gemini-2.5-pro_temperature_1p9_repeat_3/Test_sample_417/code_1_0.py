def solve_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule
    by analyzing its symmetry and counting chemically non-equivalent protons.
    """

    # --- Step-by-step analysis ---

    # 1. Central 2,4,6-trimethylbenzene core
    # Due to C3 symmetry, the three methyl groups at positions 2, 4, and 6
    # are chemically equivalent. The benzene ring has no protons of its own.
    core_methyl_signals = 1
    print(f"Signals from the three equivalent methyl groups on the central benzene ring: {core_methyl_signals}")

    # 2. The -CH2- linker group
    # Each of the three equivalent arms has a -CH2- group connecting the
    # benzene ring to the chiral indazole substituent. Because the substituent
    # is chiral, the two protons on the -CH2- group are diastereotopic
    # and therefore chemically non-equivalent.
    linker_ch2_signals = 2
    print(f"Signals from the diastereotopic protons of the -CH2- linker: {linker_ch2_signals}")

    # 3. The chiral substituent: (4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl
    # This group is a single enantiomer and has no internal plane of symmetry.
    # We must count each set of non-equivalent protons within it.

    # 3a. Protons on the indazole (fused pyrazole) ring
    # The pyrazole ring moiety has two C-H protons. Since they are in an
    # asymmetric environment, they are non-equivalent.
    indazole_ch_signals = 2
    print(f"Signals from the two non-equivalent C-H protons on the pyrazole ring: {indazole_ch_signals}")

    # 3b. Protons from the 'trimethyl' part
    # There is one methyl group at position C7.
    # There is a gem-dimethyl group at position C8. These two methyl groups
    # are diastereotopic due to the chirality of the molecule.
    # Total methyl signals = 1 (from C7) + 2 (from C8) = 3
    substituent_methyl_signals = 3
    print(f"Signals from the three non-equivalent methyl groups of the substituent: {substituent_methyl_signals}")

    # 3c. Protons on the bicyclic framework ('tetrahydro-methano' part)
    # This part lacks any symmetry. All protons not related by free rotation are unique.
    # - Bridgehead C-H proton: 1 signal
    # - Methylene group 1 (e.g., at C5): 2 diastereotopic protons -> 2 signals
    # - Methylene group 2 (e.g., at C6): 2 diastereotopic protons -> 2 signals
    aliphatic_framework_signals = 1 + 2 + 2
    print(f"Signals from the aliphatic protons on the bicyclic framework (1 CH + 2 CH2): {aliphatic_framework_signals}")

    # --- Final Calculation ---
    total_signals = (core_methyl_signals +
                     linker_ch2_signals +
                     indazole_ch_signals +
                     substituent_methyl_signals +
                     aliphatic_framework_signals)

    print("\nTotal number of expected 1H NMR signals is the sum of all unique proton environments:")
    # We use print here again so the final response only contains one code block.
    # The thought process specified each value should be in the final equation.
    print(f"Total signals = {core_methyl_signals} (core CH3) + {linker_ch2_signals} (linker CH2) + {indazole_ch_signals} (indazole CH) + {substituent_methyl_signals} (subst. CH3) + {aliphatic_framework_signals} (subst. aliphatic CH/CH2)")
    print(f"Total signals = {core_methyl_signals} + {linker_ch2_signals} + {indazole_ch_signals} + {substituent_methyl_signals} + {aliphatic_framework_signals} = {total_signals}")

solve_nmr_peaks()
<<<I>>>