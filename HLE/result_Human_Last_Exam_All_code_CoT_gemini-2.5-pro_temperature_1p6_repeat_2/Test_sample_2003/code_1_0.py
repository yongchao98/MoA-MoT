def solve_key_signature_sum():
    """
    This script derives the formula for the sum of sharps in the key signatures
    of 12 notes that have each been sharped 'n' times.
    """
    
    # The initial list of 12 musical notes.
    initial_notes = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    # Rule: Number of sharps for a natural note's key signature.
    # For F, F-major (1 flat) is converted to E-sharp major. E-major has 4 sharps,
    # so E-sharp major has 4 + 7 = 11 sharps.
    base_sharps_map = {
        'C': 0, 'D': 2, 'E': 4, 'F': 11, 'G': 1, 'A': 3, 'B': 5
    }

    # Calculate S(0), the sum of sharps for the initial list (n=0).
    s0 = 0
    for note in initial_notes:
        # The natural note is the first character (e.g., 'C' for 'C#').
        base_note_char = note[0]
        
        # The number of sharps in the note's name (e.g., 1 for 'C#').
        num_sharps_in_name = note.count('#')
        
        # The total sharps in the key is: Base Sharps + 7 * (sharps in name)
        key_signature_sharps = base_sharps_map[base_note_char] + 7 * num_sharps_in_name
        s0 += key_signature_sharps

    # For any n>0, each of the 12 notes gets an additional n sharps in its name.
    # This adds 7*n sharps to each key signature.
    # Total increase = 12 * 7 * n = 84n.
    # The coefficient for 'n' is therefore 12 * 7.
    n_coefficient = 12 * 7

    # The final formula is S(n) = S(0) + 84n.
    s0_val = int(s0)
    n_coeff_val = int(n_coefficient)
    
    print("The formula for the sum of sharps, S(n), is:")
    print(f"S(n) = {s0_val} + {n_coeff_val}*n")


solve_key_signature_sum()
<<<78 + 84*n>>>