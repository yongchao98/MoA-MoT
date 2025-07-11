def count_nmr_peaks():
    """
    Calculates the expected number of 1H NMR peaks for the given molecule based on its structure and symmetry.
    
    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    """
    
    # Step 1: Analyze the central benzene core and its methyl groups.
    # The molecule has C3 symmetry. The three methyl groups at positions 2, 4, and 6
    # of the central benzene ring are chemically equivalent.
    # They will produce one single signal in the 1H NMR spectrum.
    num_signals_core_Me = 1

    # Step 2: Analyze the -CH2- linker groups.
    # Due to the C3 symmetry, the three -CH2- linker groups are equivalent.
    # However, each linker is attached to a chiral substituent. This makes the two protons
    # within each -CH2- group diastereotopic, meaning they are chemically non-equivalent.
    # Therefore, the linker protons will give rise to two distinct signals.
    num_signals_linker_CH2 = 2

    # Step 3: Analyze one of the three equivalent chiral substituents.
    # The substituent is ((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl).
    # This is a rigid, chiral structure. We must count all unique proton environments within it.
    
    # Proton on the pyrazole ring (part of the indazole system).
    signals_pyrazole_H = 1
    
    # Protons on the rigid bicyclic (camphor-like) skeleton.
    # - There is one bridgehead proton (CH group).
    signals_skeleton_CH = 1
    # - There are two methylene groups (-CH2-). In each, the two protons are diastereotopic.
    #   The two CH2 groups themselves are in different environments. So, 2 * 2 = 4 unique protons.
    signals_skeleton_CH2s = 4
    
    # Protons on the three methyl groups attached to the skeleton.
    # Due to the overall chirality of the substituent, all three methyl groups are
    # in unique chemical environments and are non-equivalent.
    signals_skeleton_MEs = 3
    
    # Total signals from one substituent arm.
    num_signals_substituent = signals_pyrazole_H + signals_skeleton_CH + signals_skeleton_CH2s + signals_skeleton_MEs

    # Step 4: Sum the signals from all parts of the molecule.
    total_signals = num_signals_core_Me + num_signals_linker_CH2 + num_signals_substituent
    
    print("Calculation of expected 1H NMR peaks:")
    print("-" * 40)
    print(f"Signals from core methyl groups (equivalent by C3 symmetry): {num_signals_core_Me}")
    print(f"Signals from linker -CH2- groups (diastereotopic protons): {num_signals_linker_CH2}")
    print(f"Signals from one substituent arm: {num_signals_substituent}")
    print(f"  - Pyrazole C-H: {signals_pyrazole_H}")
    print(f"  - Skeleton C-H: {signals_skeleton_CH}")
    print(f"  - Skeleton -CH2- groups: {signals_skeleton_CH2s}")
    print(f"  - Skeleton -CH3 groups: {signals_skeleton_MEs}")
    print("-" * 40)
    
    # The final equation as requested.
    print(f"Final equation for total peaks:")
    print(f"{num_signals_core_Me} + {num_signals_linker_CH2} + {num_signals_substituent} = {total_signals}")
    print(f"\nThe total number of expected peaks is {total_signals}.")

count_nmr_peaks()