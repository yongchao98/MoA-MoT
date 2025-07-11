import sys

def solve_nmr_peaks():
    """
    This function calculates the expected number of 1H NMR peaks for the given molecule.
    The calculation is based on the molecule's symmetry and the chemical environment of each proton.
    """
    # Plan:
    # 1. Analyze the molecule's C3 symmetry.
    # 2. Count signals from the central trimethylbenzene core.
    # 3. Count signals from the -CH2- linker.
    # 4. Count signals from one of the chiral substituent arms.
    # 5. Sum the counts for the total number of signals.

    print("Analyzing the structure of 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene...")
    print("-" * 30)

    # Part 1: Central Mesitylene Core
    # The three methyl groups on the central benzene ring (at positions 2, 4, 6)
    # are chemically equivalent due to the molecule's C3 rotational symmetry.
    signals_core_Me = 1
    print(f"Signals from core methyl groups: {signals_core_Me}")

    # Part 2: Methylene Linkers (-CH2-)
    # The three -CH2- linkers are equivalent due to C3 symmetry. Within one linker,
    # the two protons are diastereotopic because they are adjacent to a chiral center.
    # They are chemically non-equivalent.
    signals_linker_CH2 = 2
    print(f"Signals from one -CH2- linker group: {signals_linker_CH2}")

    # Part 3: One Substituent Arm (Camphor-Indazole unit)
    # This part is chiral and has no internal symmetry. All unique protons give distinct signals.
    print("\nAnalyzing one substituent arm:")
    # H-3 proton on the pyrazole ring
    signals_H3 = 1
    # Bridgehead H-4 proton
    signals_H4 = 1
    # Diastereotopic protons of the C5-CH2 group
    signals_CH2_5 = 2
    # Diastereotopic protons of the C6-CH2 group
    signals_CH2_6 = 2
    # Methyl group at the C7 bridgehead
    signals_Me_7 = 1
    # Diastereotopic geminal methyl groups at the C8 bridge
    signals_Me_8 = 2
    
    substituent_signals = [signals_H3, signals_H4, signals_CH2_5, signals_CH2_6, signals_Me_7, signals_Me_8]
    substituent_labels = ["H-3", "H-4", "C5-H2", "C6-H2", "C7-Me", "C8-Me2"]
    
    for label, count in zip(substituent_labels, substituent_signals):
        print(f"Signals from {label}: {count}")

    # Part 4: Calculate Total
    total_signals = signals_core_Me + signals_linker_CH2 + sum(substituent_signals)
    
    print("-" * 30)
    print("Calculating the total number of expected peaks:")
    
    # Constructing the final equation string as requested
    equation_parts = [
        str(signals_core_Me),
        str(signals_linker_CH2),
        str(signals_H3),
        str(signals_H4),
        str(signals_CH2_5),
        str(signals_CH2_6),
        str(signals_Me_7),
        str(signals_Me_8)
    ]
    equation_str = " + ".join(equation_parts)
    
    print(f"Total Peaks = {equation_str} = {total_signals}")

solve_nmr_peaks()
<<<H>>>