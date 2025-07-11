def solve_nmr_peaks():
    """
    This script calculates the number of expected 1H NMR peaks for the molecule
    1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    by analyzing its structure and symmetry.
    """

    print("Step-by-step determination of the number of 1H NMR peaks:")
    print("-" * 60)

    # 1. Central Core Analysis
    print("1. Analysis of the central 2,4,6-trimethylbenzene core:")
    core_peaks = 1
    print(f"   - The three methyl groups on the central ring are equivalent due to C3 symmetry.")
    print(f"   - Contribution from core methyl groups: {core_peaks} peak.")
    print("-" * 60)

    # 2. Substituent Arm Analysis
    print("2. Analysis of one of the three identical substituent arms:")
    print("   - Due to C3 symmetry, we only need to analyze one arm.")

    # 2a. Methylene Linker
    linker_peaks = 2
    print("   - Part a: The -CH2- linker.")
    print("     - The two protons of this CH2 group are diastereotopic because they are adjacent to a large chiral center.")
    print(f"     - Contribution from the linker: {linker_peaks} peaks.")

    # 2b. Chiral Moiety
    print("   - Part b: The chiral (camphor-derived) moiety.")
    print("     - This fragment is chiral and has no internal symmetry, making non-interconvertible protons unique.")
    chiral_moiety_methyl_peaks = 3
    chiral_moiety_pyrazole_h_peaks = 1
    chiral_moiety_bridgehead_ch_peaks = 1
    chiral_moiety_framework_ch2_peaks = 4 # Two CH2 groups, each with 2 diastereotopic protons
    chiral_moiety_total_peaks = (chiral_moiety_methyl_peaks +
                                 chiral_moiety_pyrazole_h_peaks +
                                 chiral_moiety_bridgehead_ch_peaks +
                                 chiral_moiety_framework_ch2_peaks)
    print(f"       - 3 unique methyl groups -> {chiral_moiety_methyl_peaks} peaks")
    print(f"       - 1 proton on the pyrazole ring -> {chiral_moiety_pyrazole_h_peaks} peak")
    print(f"       - 1 proton at a bridgehead -> {chiral_moiety_bridgehead_ch_peaks} peak")
    print(f"       - 2 CH2 groups (4 unique protons) -> {chiral_moiety_framework_ch2_peaks} peaks")
    print(f"     - Total contribution from the chiral moiety: {chiral_moiety_total_peaks} peaks.")
    print("-" * 60)

    # 3. Final Calculation
    print("3. Final Calculation:")
    total_peaks = core_peaks + linker_peaks + chiral_moiety_total_peaks
    print("   Total Peaks = (Core Peaks) + (Linker Peaks) + (Chiral Moiety Peaks)")
    print(f"   Total Peaks = {core_peaks} + {linker_peaks} + {chiral_moiety_total_peaks} = {total_peaks}")

solve_nmr_peaks()
<<<H>>>