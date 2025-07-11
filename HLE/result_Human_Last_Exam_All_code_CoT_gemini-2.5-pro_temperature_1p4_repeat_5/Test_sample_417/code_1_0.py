def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule by analyzing its structure.
    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    """
    
    # Step 1: Count signals from the central 2,4,6-trimethylbenzene core.
    # The three methyl groups are equivalent due to C3 rotational symmetry.
    signals_core = 1
    print(f"Analysis of the central core:")
    print(f"The molecule has C3 symmetry, making the three methyl groups on the central benzene ring equivalent.")
    print(f"Number of signals from core = {signals_core}")
    print("-" * 30)

    # Step 2: Count signals from one of the three identical arms.
    # The three arms are equivalent due to C3 symmetry. We only need to analyze one.
    # The arm is: -CH2-R, where R is the chiral camphopyrazole group.
    
    print(f"Analysis of one substituent arm:")
    # The -CH2- linker's protons are diastereotopic because they are adjacent to a chiral center.
    signals_linker_ch2 = 2
    print(f"The -CH2- linker has 2 diastereotopic protons = {signals_linker_ch2} signals")
    
    # The chiral group R has several non-equivalent protons.
    # Pyrazole C-H proton:
    signals_pyrazole_h = 1
    print(f"The pyrazole ring has 1 C-H proton = {signals_pyrazole_h} signal")

    # Methyl groups on the camphor skeleton (all are non-equivalent):
    signals_methyls = 3
    print(f"The 3 methyl groups on the rigid chiral skeleton are non-equivalent = {signals_methyls} signals")

    # Protons on the saturated bicyclic frame:
    # One bridgehead C-H:
    signals_bridgehead_ch = 1
    print(f"The bridgehead has 1 C-H proton = {signals_bridgehead_ch} signal")
    
    # Two CH2 groups, with diastereotopic protons in each (2+2):
    signals_skeletal_ch2 = 4
    print(f"The skeleton has 2 CH2 groups, each with 2 diastereotopic protons (2+2) = {signals_skeletal_ch2} signals")

    signals_arm = signals_linker_ch2 + signals_pyrazole_h + signals_methyls + signals_bridgehead_ch + signals_skeletal_ch2
    print(f"\nTotal signals from one arm = {signals_linker_ch2} + {signals_pyrazole_h} + {signals_methyls} + {signals_bridgehead_ch} + {signals_skeletal_ch2} = {signals_arm}")
    print("-" * 30)

    # Step 3: Sum the signals from the core and one arm.
    total_signals = signals_core + signals_arm
    
    print("Final Calculation:")
    print(f"Total signals = (Signals from Core) + (Signals from one Arm)")
    print(f"Total signals = {signals_core} + {signals_arm} = {total_signals}")
    print("-" * 30)
    print(f"The total number of expected peaks is {total_signals}.")

count_nmr_peaks()