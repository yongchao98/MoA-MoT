def solve_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule.
    
    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene.

    The overall molecule has C3 symmetry. We count the signals from the central core
    and one of the three identical arms.
    """
    
    # 1. Signals from the central 2,4,6-trimethylbenzene core.
    # The three methyl groups are chemically equivalent due to C3 symmetry.
    core_ch3_signals = 1
    
    # 2. Signals from the linking -CH2- group.
    # This methylene group is attached to a chiral center (the entire substituent).
    # Its two protons are diastereotopic and thus chemically non-equivalent.
    linker_ch2_signals = 2
    
    # 3. Signals from one heterocyclic substituent arm.
    # This arm is chiral and has no internal symmetry plane.
    
    # 3a. Proton on the pyrazole ring at position C-3. This is a unique vinylic/aromatic proton.
    pyrazole_ch_signals = 1
    
    # 3b. Protons on the camphor-derived bicyclic system.
    
    # Methyl groups: There are three methyl groups in different environments 
    # (one bridgehead, two geminal which are syn/anti). They are all non-equivalent.
    fragment_ch3_signals = 3
    
    # Bridgehead methine proton (-CH-) at C-4. This is a single, unique proton.
    fragment_ch_signals = 1
    
    # Methylene protons (-CH2-) at C-5 and C-6.
    # The protons in each group (exo/endo) are diastereotopic.
    # The C-5 and C-6 positions are also non-equivalent.
    # This results in four distinct proton signals (H-5exo, H-5endo, H-6exo, H-6endo).
    fragment_ch2_signals = 4
    
    # 4. Sum all the signals.
    total_signals = (core_ch3_signals + 
                     linker_ch2_signals + 
                     pyrazole_ch_signals + 
                     fragment_ch3_signals + 
                     fragment_ch_signals + 
                     fragment_ch2_signals)

    print("Analysis of 1H NMR signals:")
    print("-" * 30)
    print(f"Central core (-CH3) signals: {core_ch3_signals}")
    print(f"Linking (-CH2-) group signals: {linker_ch2_signals}")
    print(f"Pyrazole ring C-H signal: {pyrazole_ch_signals}")
    print(f"Fragment methyl (-CH3) signals: {fragment_ch3_signals}")
    print(f"Fragment methine (-CH-) signal: {fragment_ch_signals}")
    print(f"Fragment methylene (-CH2-) signals: {fragment_ch2_signals}")
    print("-" * 30)
    print("Total number of expected signals is the sum of these unique environments.")
    # The final print statement shows the equation with all its components.
    print(f"Total = {core_ch3_signals} + {linker_ch2_signals} + {pyrazole_ch_signals} + {fragment_ch3_signals} + {fragment_ch_signals} + {fragment_ch2_signals}")
    print(f"Total signals = {total_signals}")

solve_nmr_peaks()