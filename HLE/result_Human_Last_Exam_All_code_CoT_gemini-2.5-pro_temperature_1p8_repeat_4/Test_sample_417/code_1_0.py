# Step 1: Define a function to calculate and explain the NMR signals.
def calculate_nmr_signals():
    """
    Calculates the number of expected 1H NMR signals for the given molecule
    by analyzing its symmetry and counting the non-equivalent protons.
    """
    print("Analyzing the molecule: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene\n")
    print("The molecule has C3 symmetry. We will count the signals for one repeating unit.\n")

    # Step 2: Count signals from the central core part of the repeating unit.
    # The three -CH3 groups on the central benzene ring are equivalent due to C3 symmetry.
    signals_central_methyls = 1
    print(f"Signals from central ring's -CH3 group: {signals_central_methyls}")

    # Step 3: Count signals from the "arm" part of the repeating unit.
    print("\n--- Analyzing one substituent arm ---")

    # The -CH2- linker group's protons are diastereotopic because they are adjacent
    # to a chiral center (the entire camphopyrazole fragment).
    signals_linker_ch2 = 2
    print(f"Signals from linker -CH2- group: {signals_linker_ch2}")

    # The camphopyrazole fragment is rigid and chiral.
    # It has one proton on the pyrazole ring.
    signals_pyrazole_ch = 1
    print(f"Signals from pyrazole ring's C-H: {signals_pyrazole_ch}")
    
    # It has a saturated bicyclic (norbornane-like) skeleton.
    # This skeleton has one unique bridgehead C-H proton.
    signals_bridgehead_ch = 1
    print(f"Signals from bicyclic skeleton's bridgehead C-H: {signals_bridgehead_ch}")
    
    # The skeleton also has two -CH2- groups. In the chiral environment, all four
    # protons are chemically distinct and non-equivalent.
    signals_skeleton_ch2 = 4
    print(f"Signals from bicyclic skeleton's two -CH2- groups: {signals_skeleton_ch2}")
    
    # The fragment has three methyl groups ('7,8,8-trimethyl'). One is on a bridgehead (C7),
    # and two are a gem-dimethyl pair on the methano bridge (C8). Due to chirality,
    # all three methyl groups are in unique environments.
    signals_fragment_methyls = 3
    print(f"Signals from the three methyl groups on the fragment: {signals_fragment_methyls}")

    # Step 4: Sum all the signals to get the total.
    total_signals = (signals_central_methyls +
                     signals_linker_ch2 +
                     signals_pyrazole_ch +
                     signals_bridgehead_ch +
                     signals_skeleton_ch2 +
                     signals_fragment_methyls)
                     
    print("\n--- Final Calculation ---")
    print("The total number of expected peaks is the sum of signals from each unique environment:")
    # Using an f-string to display the equation as requested
    final_equation = (f"{signals_central_methyls} (central CH3) + "
                      f"{signals_linker_ch2} (linker CH2) + "
                      f"{signals_pyrazole_ch} (pyrazole CH) + "
                      f"{signals_bridgehead_ch} (skeleton CH) + "
                      f"{signals_skeleton_ch2} (skeleton CH2s) + "
                      f"{signals_fragment_methyls} (fragment CH3s)")
    print(f"Total Signals = {final_equation}")
    print(f"Total Signals = {total_signals}")

# Execute the function
calculate_nmr_signals()