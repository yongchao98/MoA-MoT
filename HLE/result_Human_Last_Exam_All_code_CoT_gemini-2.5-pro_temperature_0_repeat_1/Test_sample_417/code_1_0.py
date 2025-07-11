def count_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule by analyzing its structure.

    The molecule is: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene.

    The analysis proceeds by identifying chemically non-equivalent protons based on the molecule's symmetry and structure.
    """

    # Step 1: Analyze the central 2,4,6-trimethylbenzene core.
    # The molecule has a C3 axis of symmetry. This makes the three methyl groups
    # on the central benzene ring chemically equivalent.
    # They will produce a single signal.
    signals_central_methyls = 1

    # Step 2: Analyze the methylene (-CH2-) linker.
    # The -CH2- group connects the central ring to the chiral camphor-pyrazole unit.
    # Protons on a CH2 group adjacent to a chiral center are diastereotopic,
    # meaning they are chemically non-equivalent and will produce two separate signals.
    signals_linker_ch2 = 2

    # Step 3: Analyze the pyrazole ring proton.
    # The pyrazole ring within the indazole system has a single proton on one of its carbons.
    # This proton is in a unique environment and will give one signal.
    signals_pyrazole_ch = 1

    # Step 4: Analyze the protons on the rigid, chiral camphor-derived skeleton.
    # This skeleton is asymmetric, so we count each set of non-equivalent protons.

    # Step 4a: Methyl groups on the camphor skeleton.
    # There are three methyl groups on the camphor skeleton. One is at a bridgehead,
    # and two are geminal (attached to the same carbon). Due to the chirality and
    # lack of local symmetry, all three methyl groups are in different chemical
    # environments and are thus non-equivalent.
    signals_camphor_methyls = 3

    # Step 4b: Aliphatic C-H and C-H2 protons on the camphor skeleton.
    # - There is one proton on a bridgehead carbon (-CH-). This is one signal.
    # - There are two methylene groups (-CH2-). In the rigid bicyclic system,
    #   the two protons within each CH2 group are diastereotopic (endo/exo), and the
    #   two CH2 groups are also in different environments. This results in four
    #   unique proton signals.
    # Total signals from this part = 1 (from -CH-) + 4 (from two -CH2-) = 5.
    signals_camphor_aliphatic = 5

    # Step 5: Sum all the signals to find the total number of expected peaks.
    total_signals = (signals_central_methyls +
                     signals_linker_ch2 +
                     signals_pyrazole_ch +
                     signals_camphor_methyls +
                     signals_camphor_aliphatic)

    # Print the breakdown of the calculation as requested.
    print("Calculation of expected 1H NMR peaks:")
    print(f"Central core methyls: {signals_central_methyls}")
    print(f"Methylene (-CH2-) linkers: {signals_linker_ch2}")
    print(f"Pyrazole ring C-H: {signals_pyrazole_ch}")
    print(f"Camphor skeleton methyls: {signals_camphor_methyls}")
    print(f"Camphor skeleton aliphatic C-H and C-H2: {signals_camphor_aliphatic}")
    print("-" * 20)
    # The final equation showing how the total is calculated.
    print(f"Total peaks = {signals_central_methyls} + {signals_linker_ch2} + {signals_pyrazole_ch} + {signals_camphor_methyls} + {signals_camphor_aliphatic} = {total_signals}")

# Execute the function to print the result.
count_nmr_peaks()