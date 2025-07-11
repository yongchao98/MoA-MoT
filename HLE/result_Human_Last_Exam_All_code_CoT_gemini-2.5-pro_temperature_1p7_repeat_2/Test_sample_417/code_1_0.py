def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its structure and symmetry.

    The molecule is 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene.

    Analysis Steps:
    1.  Symmetry: The molecule has C3 rotational symmetry. This means the 3 main substituent arms are equivalent,
        and the 3 methyl groups on the central benzene ring are equivalent. We analyze 1/3 of the molecule.
    2.  Central Core (2,4,6-trimethylbenzene part): The three methyl groups are equivalent due to C3 symmetry.
    3.  Linker (-CH2-): The three linker groups are equivalent. The two protons on each linker are diastereotopic
        because they are adjacent to a chiral center, making them non-equivalent.
    4.  Camphor-Indazole Moiety: We count the non-equivalent protons in one chiral unit.
    """

    # Signals from the central 2,4,6-trimethylbenzene core
    signals_core_methyls = 1

    # Signals from the -CH2- linker groups
    # The two protons are diastereotopic due to the adjacent chiral group.
    signals_linker_ch2 = 2

    # Signals from one (4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazole unit
    # This is a chiral, rigid bicyclic system.
    signals_c3_h = 1           # Proton on the pyrazole ring
    signals_c4_h = 1           # Bridgehead proton
    signals_c5_h2 = 2          # Diastereotopic protons
    signals_c6_h2 = 2          # Diastereotopic protons
    signals_c7_methyl = 1      # Methyl on a chiral bridgehead
    signals_c8_gem_dimethyl = 2 # Diastereotopic methyls adjacent to a chiral center

    # Sum of signals from the indazole moiety
    total_signals_indazole_unit = (signals_c3_h + signals_c4_h + signals_c5_h2 +
                                   signals_c6_h2 + signals_c7_methyl + signals_c8_gem_dimethyl)

    # Total number of signals for the entire molecule
    total_peaks = signals_core_methyls + signals_linker_ch2 + total_signals_indazole_unit

    print("Step-by-step calculation of 1H NMR peaks:")
    print(f"1. Signals from core methyl groups (3 equivalent CH3): {signals_core_methyls}")
    print(f"2. Signals from linker groups (3 equivalent CH2, each with 2 diastereotopic H): {signals_linker_ch2}")
    print("3. Signals from one camphor-indazole unit:")
    print(f"   - C3-H: {signals_c3_h}")
    print(f"   - C4-H: {signals_c4_h}")
    print(f"   - C5-H2: {signals_c5_h2}")
    print(f"   - C6-H2: {signals_c6_h2}")
    print(f"   - C7-CH3: {signals_c7_methyl}")
    print(f"   - C8-(CH3)2: {signals_c8_gem_dimethyl}")
    print(f"   Total for indazole unit: {total_signals_indazole_unit}")
    print("-" * 30)
    
    # The final print out showing the addition
    print("Total peaks = (Core Methyls) + (Linker CH2) + (Indazole Unit Total)")
    print(f"Total peaks = {signals_core_methyls} + {signals_linker_ch2} + {total_signals_indazole_unit}")
    print(f"Total peaks = {total_peaks}")


count_nmr_peaks()