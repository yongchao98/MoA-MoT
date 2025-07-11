def count_nmr_peaks():
    """
    Calculates the expected number of 1H NMR peaks for the given molecule by analyzing its symmetry and counting unique proton environments.
    """
    # The molecule possesses C3 rotational symmetry. We analyze the protons of the central ring
    # and one of the three identical substituent arms.

    # 1. Protons on the central 2,4,6-trimethylbenzene core:
    # The three methyl groups are equivalent due to C3 symmetry.
    core_methyl_peaks = 1

    # 2. Protons in one substituent arm:
    # The arm is -CH2-[chiral camphor-indazole moiety].

    # The benzylic -CH2- group is adjacent to a chiral center, making its two protons diastereotopic and non-equivalent.
    benzylic_ch2_peaks = 2

    # The chiral camphor-indazole moiety is rigid and asymmetric.
    # a. Proton on the pyrazole ring (C3-H): Unique environment.
    pyrazole_ch_peak = 1

    # b. Three methyl groups on the camphor skeleton:
    #    - Two are geminal and are diastereotopic.
    #    - One is on a bridgehead and is unique.
    #    All three are in different chemical environments.
    camphor_methyl_peaks = 3

    # c. Bridgehead proton (CH): Unique environment.
    bridgehead_ch_peak = 1

    # d. Two CH2 groups in the saturated bicyclic ring:
    #    In this rigid, chiral system, all four protons are non-equivalent.
    ring_ch2_peaks = 4

    # 3. Sum all unique proton environments.
    total_peaks = (core_methyl_peaks +
                   benzylic_ch2_peaks +
                   pyrazole_ch_peak +
                   camphor_methyl_peaks +
                   bridgehead_ch_peak +
                   ring_ch2_peaks)

    # Output the breakdown of the calculation.
    print("Step-by-step calculation of 1H NMR peaks:")
    print(f"1. Central core methyl groups (equivalent by C3 symmetry): {core_methyl_peaks} peak")
    print(f"2. Benzylic -CH2- group (diastereotopic protons): {benzylic_ch2_peaks} peaks")
    print(f"3. Pyrazole ring C-H proton: {pyrazole_ch_peak} peak")
    print(f"4. Camphor skeleton methyl groups (all non-equivalent): {camphor_methyl_peaks} peaks")
    print(f"5. Camphor skeleton bridgehead C-H proton: {bridgehead_ch_peak} peak")
    print(f"6. Camphor skeleton ring -CH2- protons (all non-equivalent): {ring_ch2_peaks} peaks")
    print("-" * 30)

    # Print the final equation with each number.
    print("Final Equation:")
    print(f"Total Peaks = {core_methyl_peaks} + {benzylic_ch2_peaks} + {pyrazole_ch_peak} + {camphor_methyl_peaks} + {bridgehead_ch_peak} + {ring_ch2_peaks}")
    print(f"\nTotal number of expected peaks: {total_peaks}")

count_nmr_peaks()