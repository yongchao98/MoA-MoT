def derive_music_theory_formula():
    """
    This program derives the formula S(n) for the total sum of sharps based on the problem description.
    It breaks the problem down into sharps from tonic names and sharps from key signatures,
    calculates the constant components, and assembles the final formula.
    """

    print("Deriving the formula S(n) for the total sum of sharps.")
    print("The total sum S(n) is the sum of two components calculated across 12 keys:\n")
    print("1. The total number of sharps in the names of the tonic notes.")
    print("2. The total number of sharps in the corresponding key signatures.\n")

    # Part 1: Calculating sharps from the names of the tonics.
    # The initial 12 notes are C, C#, D, D#, E, F, F#, G, G#, A, A#, B.
    # The notes with one sharp in their name are C#, D#, F#, G#, A#.
    initial_sharps_in_names = 5
    print("--- Part 1: Sharps from Tonic Names ---")
    print(f"The sum of sharps in the initial 12 note names (C#, D#, etc.) is a constant: {initial_sharps_in_names}.")
    print("Each of the 12 notes is then sharped 'n' times, adding 'n' sharps to each name.")
    print("This adds a variable component of 12 * n to the total sum.")
    print(f"Total sharps from all tonic names = {initial_sharps_in_names} + 12 * n\n")

    # Part 2: Calculating sharps from the key signatures.
    # The number of sharps for a tonic with pitch class 'p' (C=0) is (p * 7) % 12.
    # Summing this from p=0 to 11 gives a constant value, as the set of pitch
    # classes is just rotated by 'n'.
    sharps_in_signatures_sum = sum((p * 7) % 12 for p in range(12))
    print("--- Part 2: Sharps from Key Signatures ---")
    print("The sum of sharps in the 12 major key signatures is also a constant, independent of 'n'.")
    print(f"This is because the set of 12 tonic pitch classes is always {0, 1, ..., 11}, just rotated.")
    print(f"The sum is calculated as Sum_{{p=0 to 11}} (p * 7) % 12, which equals {sharps_in_signatures_sum}.\n")
    
    # Combine the parts to get the final formula S(n).
    final_constant_term = initial_sharps_in_names + sharps_in_signatures_sum
    print("--- Final Formula ---")
    print("Combining the components gives the total sum S(n):")
    print("S(n) = (Sharps from names) + (Sharps from signatures)")
    print(f"S(n) = ({initial_sharps_in_names} + 12 * n) + {sharps_in_signatures_sum}")
    print(f"S(n) = 12 * n + ({initial_sharps_in_names} + {sharps_in_signatures_sum})")
    
    print("\nThe final simplified formula, with each number shown, is:")
    print(f"S(n) = 12 * n + {final_constant_term}")


# Execute the derivation
derive_music_theory_formula()