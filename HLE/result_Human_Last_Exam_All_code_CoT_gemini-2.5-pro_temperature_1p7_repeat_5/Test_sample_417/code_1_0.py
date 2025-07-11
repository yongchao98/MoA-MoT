def count_nmr_signals():
    """
    This function calculates the number of expected 1H NMR signals for
    1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    based on its molecular structure and symmetry.
    """

    # Due to C3 symmetry, we analyze the central core and one of the three identical arms.

    # 1. Signals from the central 2,4,6-trimethylbenzene core.
    # The three methyl groups are equivalent.
    core_signals = 1

    # 2. Signals from one substituent arm.
    # The arm consists of a -CH2- linker and a chiral indazole unit.
    
    # 2a. The -CH2- linker has two diastereotopic protons.
    linker_signals = 2

    # 2b. The chiral indazole moiety has a fixed, asymmetric structure.
    # Counting its non-equivalent protons:
    # - 1 vinylic H on the pyrazole ring
    # - 1 bridgehead H
    # - 2 diastereotopic Hs on one CH2 group
    # - 2 diastereotopic Hs on another CH2 group
    # - 1 methyl group
    # - 2 diastereotopic methyls on the gem-dimethyl bridge
    indazole_signals = 1 + 1 + 2 + 2 + 1 + 2

    # 3. Summing up all the signals.
    total_signals = core_signals + linker_signals + indazole_signals
    
    print("The total number of expected 1H NMR signals is calculated by summing the signals from each unique part of the molecule:")
    print("Number of signals from core methyls: 1")
    print("Number of signals from the CH2 linker: 2")
    print("Number of signals from the chiral indazole moiety: 9")
    print("\nFinal calculation:")
    print(f"{core_signals} + {linker_signals} + {indazole_signals} = {total_signals}")

count_nmr_signals()