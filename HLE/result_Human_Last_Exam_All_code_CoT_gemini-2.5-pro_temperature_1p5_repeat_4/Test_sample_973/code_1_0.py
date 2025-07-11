import re

def solve_musical_puzzle():
    """
    This script determines the final sung note of "Happy Birthday" based on
    the provided chord progression.
    """

    # The full chord progression provided for the song.
    progression = [
        "Cm7 F7(9)",
        "Bm7 E7(9)",
        "Am7 D7(9)",
        "Abm7 Db7(9)",
        "Ebm7 Ab7(9)",
        "Bm7 E7(9)",
        "Cm7 F7(9)"
    ]

    # A mapping of musical notes to numeric values (semitones from C).
    note_to_value = {
        'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3, 'E': 4, 'F': 5,
        'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11
    }
    
    # Mapping of values back to note names (preferring flats for common jazz keys).
    value_to_note = {
        0: 'C', 1: 'Db', 2: 'D', 3: 'Eb', 4: 'E', 5: 'F', 6: 'Gb', 7: 'G', 8: 'Ab', 9: 'A', 10: 'Bb', 11: 'B'
    }

    print("Step 1: Identifying the final chords from the progression.")
    # The song structure implies the final chords provided are for the final "happy birthday".
    # The final chord pair for "bir-day" is the last one in the list.
    final_pair_str = progression[-1]
    final_ii_chord, final_V_chord = final_pair_str.split()
    print(f"The last phrase uses '{final_ii_chord}' on 'bir...' and '{final_V_chord}' on '...day'.")
    print("-" * 20)

    print("Step 2: Analyzing the harmony.")
    print(f"The chord pair {final_ii_chord} - {final_V_chord} is a 'ii-V' progression.")
    print("This progression creates a strong musical expectation to resolve to a tonic ('I') key.")
    print("-" * 20)
    
    print("Step 3: Calculating the tonic key.")
    # Extract the root note from the V chord, e.g., 'F' from 'F7(9)'.
    v_chord_root = re.match(r"[A-G][b#]?", final_V_chord).group(0)
    v_chord_value = note_to_value[v_chord_root]

    print("The tonic key is a perfect fifth (7 semitones) below the root of the V chord.")
    print(f"The equation to find the tonic's value is: (Value of '{v_chord_root}' - 7) % 12")
    
    # Calculate the tonic value
    tonic_value = (v_chord_value - 7 + 12) % 12
    
    print(f"The calculation is: ({v_chord_value} - 7) % 12 = {tonic_value}")
    
    tonic_note = value_to_note[tonic_value]
    print(f"The note with value {tonic_value} is '{tonic_note}'. So, the key is {tonic_note} Major.")
    print("-" * 20)

    print("Step 4: Determining the final melody note.")
    print("The standard melody for 'Happy Birthday to You' always ends on the tonic note of the key.")
    print(f"Since the key is {tonic_note} Major, the final note sung on the word 'you' is the tonic.")
    print("\nTherefore, the concluding note is:")
    
    final_answer = tonic_note
    print(final_answer)

solve_musical_puzzle()
<<<Bb>>>