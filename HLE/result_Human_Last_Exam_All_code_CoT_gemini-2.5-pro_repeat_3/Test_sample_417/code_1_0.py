def solve_nmr_peaks():
    """
    This function analyzes the structure of 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-
    4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    to determine the number of expected 1H NMR signals.
    """

    print("Step-by-step analysis of the number of 1H NMR peaks:")
    print("-" * 50)

    # Step 1: Analyze the central core and linkers
    print("1. Analysis of the Central Core and Linker Groups:")
    print("   - The molecule possesses a C3 rotational symmetry axis, making the three arms and three core methyl groups equivalent.")
    
    # Core methyl groups
    core_methyl_signals = 1
    print(f"   - The three equivalent methyl groups (-CH3) on the benzene ring give {core_methyl_signals} signal.")
    
    # Linker methylene groups
    linker_ch2_signals = 2
    print(f"   - The protons on each of the three equivalent methylene (-CH2-) linkers are diastereotopic due to the adjacent chiral substituent.")
    print(f"   - This gives {linker_ch2_signals} signals for the linker protons.")

    total_core_signals = core_methyl_signals + linker_ch2_signals
    print(f"   - Total signals from the core and linkers = {core_methyl_signals} + {linker_ch2_signals} = {total_core_signals}")
    print("-" * 50)

    # Step 2: Analyze one of the substituent arms
    print("2. Analysis of one Substituent Arm:")
    print("   - The arm is a rigid, chiral structure with no internal symmetry.")
    
    # Protons on the arm
    pyrazole_h_signals = 1
    print(f"   - One proton on the pyrazole ring gives {pyrazole_h_signals} signal.")
    
    methyl_signals = 3
    print(f"   - Three non-equivalent methyl groups in the bicyclic system give {methyl_signals} signals.")

    bridgehead_h_signals = 1
    print(f"   - One bridgehead proton gives {bridgehead_h_signals} signal.")

    aliphatic_ch2_signals = 4
    print(f"   - Four non-equivalent protons from the two aliphatic CH2 groups give {aliphatic_ch2_signals} signals.")

    total_arm_signals = pyrazole_h_signals + methyl_signals + bridgehead_h_signals + aliphatic_ch2_signals
    print(f"   - Total signals from one arm = {pyrazole_h_signals} + {methyl_signals} + {bridgehead_h_signals} + {aliphatic_ch2_signals} = {total_arm_signals}")
    print("-" * 50)

    # Step 3: Calculate the total number of signals
    print("3. Calculation of Total Number of Signals:")
    total_signals = total_core_signals + total_arm_signals
    print("   - Total = (Signals from Core/Linkers) + (Signals from one Arm)")
    
    # Printing each number in the final equation
    print(f"   - Final Calculation: {total_core_signals} + {total_arm_signals} = {total_signals}")
    
    print("-" * 50)
    print(f"The total number of expected 1H NMR peaks is {total_signals}.")

solve_nmr_peaks()