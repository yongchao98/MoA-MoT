def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its structure and symmetry.

    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    """
    
    # Analysis of the molecule's structure and symmetry:
    # The molecule has a C3 axis of symmetry, making the three large substituents equivalent.
    # We count the unique proton signals from the central core and one of the substituents.

    # 1. Signals from the central 2,4,6-trimethylbenzene core and linkers
    
    # The three methyl groups on the benzene ring (at C2, C4, C6) are equivalent by symmetry.
    mesitylene_me_signals = 1
    
    # The three -CH2- linker groups (at C1, C3, C5) are equivalent by symmetry.
    # However, each CH2 is attached to a chiral center, making its two protons diastereotopic (non-equivalent).
    linker_ch2_signals = 2
    
    # 2. Signals from one chiral substituent unit
    # This unit is ((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)
    # It is a rigid, chiral structure, so most protons are in unique environments.
    
    # There are three non-equivalent methyl groups in the camphor-derived skeleton.
    substituent_me_signals = 3
    
    # There is one proton on the pyrazole ring.
    substituent_pyrazole_h_signals = 1
    
    # There is one methine proton at a bridgehead of the bicyclic system.
    substituent_bridgehead_h_signals = 1
    
    # There are two CH2 groups in the bicyclic system. All four protons are diastereotopic.
    substituent_ring_ch2_signals = 4
    
    # 3. Summing all the signals
    total_signals = (mesitylene_me_signals + 
                     linker_ch2_signals + 
                     substituent_me_signals + 
                     substituent_pyrazole_h_signals +
                     substituent_bridgehead_h_signals +
                     substituent_ring_ch2_signals)

    print("Step-by-step calculation of 1H NMR peaks:")
    print("-" * 40)
    print(f"Signals from equivalent mesitylene methyl groups: {mesitylene_me_signals}")
    print(f"Signals from equivalent but diastereotopic linker -CH2- groups: {linker_ch2_signals}")
    print(f"Signals from non-equivalent substituent methyl groups: {substituent_me_signals}")
    print(f"Signals from substituent pyrazole ring C-H: {substituent_pyrazole_h_signals}")
    print(f"Signals from substituent bridgehead C-H: {substituent_bridgehead_h_signals}")
    print(f"Signals from substituent diastereotopic ring -CH2- protons: {substituent_ring_ch2_signals}")
    print("-" * 40)
    
    # To meet the output requirement, we print the equation explicitly with the numbers.
    print(f"Total number of peaks = {mesitylene_me_signals} + {linker_ch2_signals} + {substituent_me_signals} + {substituent_pyrazole_h_signals} + {substituent_bridgehead_h_signals} + {substituent_ring_ch2_signals} = {total_signals}")

count_nmr_peaks()