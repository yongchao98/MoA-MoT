def calculate_nmr_peaks():
    """
    Calculates the expected number of 1H NMR peaks for the given molecule based on its symmetry and structure.
    """

    # Step 1: Analyze the central core (2,4,6-trimethylbenzene).
    # The three methyl groups are equivalent due to the molecule's C3 symmetry.
    # The other ring positions are substituted, so there are no aromatic protons.
    signals_from_core = 1
    print(f"Number of signals from the central core's methyl groups: {signals_from_core}")

    # Step 2: Analyze one of the three identical substituent arms.
    # The arm is chiral, so we must consider diastereotopicity.

    # 2a. The benzylic methylene linker (-CH2-):
    # Its two protons are diastereotopic because they are adjacent to a chiral system.
    signals_from_linker = 2
    print(f"Number of signals from the linker -CH2- group: {signals_from_linker}")

    # 2b. The proton on the pyrazole ring (part of the indazole system):
    # There is one unique proton on the five-membered pyrazole ring.
    signals_from_pyrazole_H = 1
    print(f"Number of signals from the pyrazole ring proton: {signals_from_pyrazole_H}")

    # 2c. The methyl groups on the rigid bicyclic frame:
    # - One unique methyl group at position C7.
    # - Two diastereotopic methyl groups at position C8.
    signals_from_camphor_Me = 1 + 2
    print(f"Number of signals from the three methyl groups on the bicyclic frame: {signals_from_camphor_Me}")

    # 2d. The CH and CH2 protons on the rigid bicyclic frame:
    # - One unique bridgehead proton at C4.
    # - Two diastereotopic protons on the C5 methylene group.
    # - Two diastereotopic protons on the C6 methylene group.
    signals_from_camphor_frame_protons = 1 + 2 + 2
    print(f"Number of signals from the CH/CH2 protons on the bicyclic frame: {signals_from_camphor_frame_protons}")

    # Step 3: Sum the signals from one arm.
    signals_from_one_arm = signals_from_linker + signals_from_pyrazole_H + signals_from_camphor_Me + signals_from_camphor_frame_protons
    print(f"\nTotal signals from one substituent arm: {signals_from_one_arm}")

    # Step 4: Calculate the total number of signals for the entire molecule.
    total_signals = signals_from_core + signals_from_one_arm
    
    print("\nTotal expected 1H NMR peaks = (signals from core) + (signals from one arm)")
    # Final equation as requested, showing each number
    print(f"{signals_from_core} + {signals_from_one_arm} = {total_signals}")

calculate_nmr_peaks()