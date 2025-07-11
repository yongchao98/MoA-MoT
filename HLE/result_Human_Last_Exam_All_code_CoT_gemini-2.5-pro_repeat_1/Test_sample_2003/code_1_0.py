def solve_music_formula():
    """
    This script derives the formula for the sum of sharps in the major key signatures
    of 12 notes after they have been sharpened n times (for n > 0).
    """

    # Step 1: Define the number of sharps for the 7 natural notes based on the Circle of Fifths.
    # F has -1 sharps (i.e., 1 flat).
    natural_note_sharps = {
        'C': 0, 'G': 1, 'D': 2, 'A': 3, 'E': 4, 'B': 5, 'F': -1
    }

    # Step 2: Define the 12 initial notes as a combination of a natural note and a number of sharps.
    # For example, C# is (base_note='C', initial_sharps=1).
    initial_notes = [
        ('C', 0),  # C
        ('C', 1),  # C#
        ('D', 0),  # D
        ('D', 1),  # D#
        ('E', 0),  # E
        ('F', 0),  # F
        ('F', 1),  # F#
        ('G', 0),  # G
        ('G', 1),  # G#
        ('A', 0),  # A
        ('A', 1),  # A#
        ('B', 0)   # B
    ]

    n_coefficient = 0
    constant_term = 0

    # Step 3: For each note, calculate its contribution to the final formula S(n) = a*n + b.
    # A note X with k sharps, sharpened n times, becomes X with (k+n) sharps.
    # The number of sharps in its key signature is NumSharps(X) + 7*(k+n).
    # This can be expanded to 7n + (NumSharps(X) + 7k).
    for base_note, initial_sharps in initial_notes:
        # Each note adds 7 to the coefficient of n.
        n_coefficient += 7
        
        # Calculate the constant part for this note.
        base_sharps = natural_note_sharps[base_note]
        constant_part = base_sharps + 7 * initial_sharps
        constant_term += constant_part

    # Step 4: Print the final derived formula.
    print("The derived formula for the sum of sharps for n > 0 is S(n) = a*n + b")
    print(f"Calculated coefficient for n (a): {n_coefficient}")
    print(f"Calculated constant term (b): {constant_term}")
    print("Final Formula:")
    # The requirement is to output each number in the final equation.
    print(f"S(n) = {n_coefficient}*n + {constant_term}")


solve_music_formula()
<<<84*n + 54>>>