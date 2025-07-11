def count_nmr_peaks():
    """
    This function calculates the expected number of 1H NMR peaks for the molecule
    1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    by analyzing its symmetry and structure.
    """
    
    # 1. Central Core Signals (2,4,6-trimethylbenzene)
    # The three methyl groups on the central benzene ring are equivalent due to C3 symmetry.
    signals_mesityl_core = 1

    # 2. Linker Signals (-CH2- group)
    # The molecule has three equivalent -CH2- linkers. Within one linker, the two protons
    # are adjacent to a chiral center (the indazole moiety), making them diastereotopic
    # and thus chemically non-equivalent.
    signals_linker = 2

    # 3. Indazole Moiety Signals
    # We analyze one of the three identical chiral indazole units.
    
    # Proton on the pyrazole ring (C3-H). It's a unique aromatic-like proton.
    signals_pyrazole_H = 1
    
    # Methyl groups on the indazole scaffold (7,8,8-trimethyl).
    # The methyl at C7 is unique.
    signals_C7_Me = 1
    # The two geminal methyls at C8 are diastereotopic due to the overall chirality.
    signals_C8_diMe = 2
    
    # Protons on the rigid bicyclic frame.
    # The single proton at the C4 bridgehead.
    signals_C4_H = 1
    # The two protons at C5 (CH2) are diastereotopic (exo and endo).
    signals_C5_H2 = 2
    # The two protons at C6 (CH2) are also diastereotopic (exo and endo).
    signals_C6_H2 = 2
    
    # 4. Total Calculation
    total_signals = (signals_mesityl_core +
                     signals_linker +
                     signals_pyrazole_H +
                     signals_C7_Me +
                     signals_C8_diMe +
                     signals_C4_H +
                     signals_C5_H2 +
                     signals_C6_H2)
    
    # Printing the breakdown of the calculation as an equation
    print("The total number of signals is calculated by summing the signals from each unique part:")
    print(f"Total Signals = {signals_mesityl_core} (central core methyls) + "
          f"{signals_linker} (linker CH2) + "
          f"{signals_pyrazole_H} (pyrazole H) + "
          f"{signals_C7_Me} (C7-Me) + {signals_C8_diMe} (C8-diMe) + "
          f"{signals_C4_H} (C4-H) + {signals_C5_H2} (C5-CH2) + {signals_C6_H2} (C6-CH2)")
          
    print(f"\nFinal Equation:")
    print(f"{signals_mesityl_core} + {signals_linker} + {signals_pyrazole_H} + {signals_C7_Me} + {signals_C8_diMe} + {signals_C4_H} + {signals_C5_H2} + {signals_C6_H2} = {total_signals}")

    print(f"\nTotal number of expected peaks in the 1H NMR spectrum: {total_signals}")

count_nmr_peaks()