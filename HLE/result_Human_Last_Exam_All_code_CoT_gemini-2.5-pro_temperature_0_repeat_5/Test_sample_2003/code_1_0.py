import sys

def solve_music_theory_sharps():
    """
    This script derives and verifies a formula for the sum of sharps of 12 major key signatures.
    """

    # Step 1: Define the musical data
    # The 12 initial notes are defined by their letter name and number of accidentals (sharps).
    initial_notes = [
        ('C', 0), ('C', 1),  # C, C#
        ('D', 0), ('D', 1),  # D, D#
        ('E', 0),            # E
        ('F', 0), ('F', 1),  # F, F#
        ('G', 0), ('G', 1),  # G, G#
        ('A', 0), ('A', 1),  # A, A#
        ('B', 0)             # B
    ]

    # Base sharps for natural major keys (C, D, E, F, G, A, B). F has -1 sharps (1 flat).
    s_base = {'C': 0, 'D': 2, 'E': 4, 'F': -1, 'G': 1, 'A': 3, 'B': 5}
    
    # Pitch class values for enharmonic conversion (C=0, C#=1, etc.).
    pitch_class_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}

    def calculate_total_sharps(n):
        """
        Calculates the total number of sharps for a given n based on the problem's rules.
        """
        total_sharps = 0
        for letter, initial_acc in initial_notes:
            # For each note in the initial list, add n sharps.
            current_acc = initial_acc + n
            
            # Calculate the "raw" number of sharps.
            # Formula: S(Note) = S_base(letter) + 7 * number_of_accidentals
            raw_sharps = s_base[letter] + 7 * current_acc
            
            # Check if the key has flats (raw_sharps < 0).
            # If so, apply the enharmonic conversion to a sharp key.
            if raw_sharps < 0:
                # The number of sharps in the enharmonic key is (p * 7) % 12,
                # where p is the note's pitch class (0-11).
                p_class = (pitch_class_map[letter] + current_acc) % 12
                final_sharps = (p_class * 7) % 12
            else:
                # If the key has 0 or more sharps, no conversion is needed.
                final_sharps = raw_sharps
            
            total_sharps += final_sharps
            
        return total_sharps

    # Step 2 & 3: Derive the formula by analyzing n=0 and n>0 cases.
    # The logic shows that the flat-to-sharp conversion only happens for the note F when n=0.
    # For n=0, S_raw(F) = -1, which is converted to 11 sharps. The difference is 12.
    # The sum of raw sharps for n=0 is 54. With the conversion, it becomes 54 - (-1) + 11 = 66.
    # For n>0, no key has flats, so the sum is the sum of raw sharps.
    # Sum S(n) = Sum(S_raw(note_i) + 7n) = (Sum S_raw(note_i)) + 12*7n = 54 + 84n.

    # Step 4: Print the final derived formula.
    print("--- Derivation and Formula ---")
    print("Let S(n) be the total number of sharps.")
    print("The formula is derived by considering two cases based on the input 'n'.\n")
    print("Case 1: n = 0")
    print("The key of F major has 1 flat. As per the rules, it is converted to its")
    print("enharmonic equivalent, E# major, which has 11 sharps.")
    print("The sum for all 12 keys is calculated to be 66.")
    
    print("\nCase 2: n > 0")
    print("When n is greater than 0, all 12 resulting keys have 0 or more sharps.")
    print("No flat-to-sharp conversion is needed. The sum follows a linear formula.")

    print("\n--- Final Simplified Formula ---")
    val_n0 = 66
    coeff_n = 84
    const_n_gt_0 = 54
    
    print(f"For n = 0: S(0) = {val_n0}")
    print(f"For n > 0: S(n) = {coeff_n}n + {const_n_gt_0}")
    print(f"The formula for n > 0 can also be simplified to: S(n) = 6 * ({coeff_n//6}n + {const_n_gt_0//6})")

    # Step 5: Verify the formula with the calculation function.
    print("\n--- Verification ---")
    for n_val in range(4):
        calculated_sum = calculate_total_sharps(n_val)
        if n_val == 0:
            formula_sum = val_n0
        else:
            formula_sum = coeff_n * n_val + const_n_gt_0
        
        print(f"For n={n_val}: Calculated Sum = {calculated_sum}, Formula Sum = {formula_sum}. Match: {calculated_sum == formula_sum}")

solve_music_theory_sharps()
# The triple angle bracket answer format is required by the system.
# It is not part of the python code's output.
sys.stdout.flush()
print("<<<S(n) = 66 for n=0; S(n) = 84n + 54 for n>0>>>")