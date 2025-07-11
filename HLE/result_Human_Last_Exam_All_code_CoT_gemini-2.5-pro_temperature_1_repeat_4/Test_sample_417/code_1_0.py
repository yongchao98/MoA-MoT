def calculate_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule by analyzing its symmetry and structure.
    """
    print("Analyzing the structure: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene")
    print("The molecule has C3 symmetry, so we count unique protons in one of the three identical repeating units.\n")

    # Part 1: Protons on the central core and linkers
    # The three methyl groups on the benzene ring are equivalent.
    signals_ar_ch3 = 1
    print(f"Signals from the three equivalent methyl groups on the central benzene ring: {signals_ar_ch3}")

    # The -CH2- linkers are equivalent, but the two protons on each linker are diastereotopic.
    signals_linker_ch2 = 2
    print(f"Signals from the diastereotopic protons of the -CH2- linkers: {signals_linker_ch2}")

    # Part 2: Protons on one substituent arm
    # The arm is a chiral, rigid camphor-indazole derivative.
    print("\nAnalyzing one substituent arm:")

    # Proton on the pyrazole ring (H-3)
    signals_pyrazole_h = 1
    print(f"- Signal from the proton on the pyrazole ring: {signals_pyrazole_h}")

    # Protons on the bicyclic framework
    signals_bridgehead_h = 1
    print(f"- Signal from the bridgehead methine proton (C4-H): {signals_bridgehead_h}")
    signals_c5_h2 = 2
    print(f"- Signals from the diastereotopic C5-H protons (endo/exo): {signals_c5_h2}")
    signals_c6_h2 = 2
    print(f"- Signals from the diastereotopic C6-H protons (endo/exo): {signals_c6_h2}")

    # Protons on the methyl groups of the arm
    signals_c7_me = 1
    print(f"- Signal from the C7-methyl group: {signals_c7_me}")
    signals_c8_me2 = 2
    print(f"- Signals from the two diastereotopic C8-methyl groups: {signals_c8_me2}")

    # Part 3: Sum all signals for the final equation and answer
    total_signals = (signals_ar_ch3 + signals_linker_ch2 + signals_pyrazole_h +
                     signals_bridgehead_h + signals_c5_h2 + signals_c6_h2 +
                     signals_c7_me + signals_c8_me2)

    print("\n---\nFinal Calculation:")
    print("The total number of peaks is the sum of all chemically distinct proton signals.")
    print(f"{signals_ar_ch3} (Ar-CH3) + {signals_linker_ch2} (-CH2-) + {signals_pyrazole_h} (H-3) + {signals_bridgehead_h} (C4-H) + {signals_c5_h2} (C5-H2) + {signals_c6_h2} (C6-H2) + {signals_c7_me} (C7-Me) + {signals_c8_me2} (C8-Me2) = {total_signals}")

    print(f"\nTotal number of expected 1H NMR peaks: {total_signals}")

calculate_nmr_peaks()