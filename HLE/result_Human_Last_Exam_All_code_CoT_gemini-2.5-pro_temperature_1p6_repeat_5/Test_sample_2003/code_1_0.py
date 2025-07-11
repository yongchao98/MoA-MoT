import sys

def solve_music_theory_problem():
    """
    This function derives a formula for the sum of sharps in the key signatures
    of 12 musical notes after they have all been sharpened 'n' times.
    """

    # Step 1: Define the 12 notes and their pitch classes (0-11).
    note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    initial_pitch_classes = list(range(12))

    # Step 2: Define the function to calculate sharps for a major key.
    # For a key starting on a note with pitch class p, the number of sharps
    # (with flat keys converted) is (7 * p) mod 12.
    def get_number_of_sharps(pitch_class):
        return (7 * pitch_class) % 12

    # Step 3: Calculate the sum for the base case (n=0).
    print("--- Derivation of the Formula ---")
    print("\nStep 1: Calculate the sum for the initial 12 notes (the case where n=0).")
    print("The 12 notes are C, C#, D, D#, E, F, F#, G, G#, A, A#, B.")
    
    sharps_per_note = []
    for p in initial_pitch_classes:
        sharps = get_number_of_sharps(p)
        sharps_per_note.append(sharps)
    
    total_sharps = sum(sharps_per_note)

    # Print the full equation for the sum
    sum_equation_str = " + ".join(map(str, sharps_per_note))
    print(f"\nThe number of sharps for the major key of each note are:")
    print(f"[{', '.join(map(str, sharps_per_note))}]")
    print(f"\nThe sum is the following equation:")
    print(f"Sum = {sum_equation_str} = {total_sharps}")

    # Step 4: Analyze the effect of 'n'.
    print("\nStep 2: Analyze the effect of sharpening all notes 'n' times.")
    print("When each of the 12 unique notes is sharpened 'n' times, the resulting set of 12 notes is simply a permutation of the original set.")
    print("For example, if n=1, C becomes C#, C# becomes D, ..., and B becomes B# (enharmonically C). The set of notes is still {C, C#, ..., B}.")
    print("Because the collection of notes is identical for any 'n', the sum of the sharps of their key signatures must also be constant.")

    # Step 5: State the final formula.
    print("\nStep 3: State the final formula in terms of n.")
    print("The sum does not depend on 'n' and is always equal to the sum calculated in Step 1.")
    
    # To satisfy the prompt of "output each number in the final equation"
    # we can express the formula S(n) as 0*n + 66.
    coefficient_n = 0
    constant_term = total_sharps
    print(f"\nThe unsimplified formula for the sum S(n) is:")
    print(f"S(n) = {coefficient_n} * n + {constant_term}")

    print("\nSimplifying this gives the final formula:")
    print(f"S(n) = {constant_term}")

solve_music_theory_problem()