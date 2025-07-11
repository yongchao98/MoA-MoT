def solve_music_formula():
    """
    This script derives the formula for the sum of sharps in key signatures.

    The logic is as follows:
    1. The number of sharps in a major key is determined by its root note's
       letter name and its accidentals. We can define a base number of sharps
       for each natural note (C,D,E,F,G,A,B), where F has -1.
    2. Adding a sharp to a tonic adds 7 sharps to its key signature.
    3. We start with a list of 12 notes: C, C#, D, ..., B.
    4. For a general case n > 0, each of these 12 notes is sharped 'n' times.
    5. The new number of sharps for each note's key signature will be its
       original number of sharps plus 7*n.
    6. The total sum is the sum of the original sharps for all 12 notes,
       plus the sum of 7*n for all 12 notes.
       Total(n) = Sum(Original_Sharps) + 12 * (7*n)
    7. We calculate the Sum(Original_Sharps) and the coefficient for n.
    """
    
    # Base sharps for natural notes (C=0, G=1, D=2, A=3, E=4, B=5, F=-1)
    base_sharps = {"C": 0, "D": 2, "E": 4, "F": -1, "G": 1, "A": 3, "B": 5}
    
    # The 12 notes from the problem description
    initial_notes = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    
    constant_term = 0
    
    # Calculate the sum of sharps for the initial 12 notes (the constant term)
    # This calculation is for the general case n>0, where the note 'F' is treated as F, not E#.
    for note in initial_notes:
        base_note = note[0]
        num_sharps_on_name = note.count('#')
        
        key_signature_sharps = base_sharps[base_note] + 7 * num_sharps_on_name
        constant_term += key_signature_sharps
        
    # For each of the 12 notes, adding 'n' sharps adds 7*n to the key signature.
    # The total contribution from 'n' is therefore 12 * 7 * n.
    n_coefficient = 12 * 7
    
    print("Deriving the formula S(n) = a*n + b for n > 0...")
    print(f"The coefficient for n, 'a', is 12 notes * 7 sharps per added sharp = {n_coefficient}")
    print(f"The constant term, 'b', is the sum of sharps for the base list (C, C#, ... B) = {constant_term}")
    print("\nThe final formula for the sum of sharps S(n) is:")
    print(f"S(n) = {n_coefficient} * n + {constant_term}")

solve_music_formula()