import re

def solve_final_note():
    """
    Calculates the final note of the song based on the provided chord progression.
    """
    # The chromatic scale. We'll use flats for consistency.
    notes = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B']

    # The chord progression from the song arrangement.
    progression = [
        "Cm7 F7(9)",
        "Bm7 E7(9)",
        "Am7 D7(9)",
        "Abm7 Db7(9)",
        "Ebm7 Ab7(9)",
        "Bm7 E7(9)",
        "Cm7 F7(9)"
    ]

    # 1. Get the final chord pair from the progression.
    final_pair = progression[-1]
    print(f"The final chord progression is: {final_pair}")

    # 2. This is a ii-V progression. The second chord is the 'V' chord.
    # We extract the root of the V chord.
    # For "Cm7 F7(9)", the V chord is F7(9) and its root is F.
    v_chord_root_str = final_pair.split(" ")[1]
    # Use regex to get the note name (handling flats and sharps if necessary)
    root_note_of_v = re.match(r'[A-G][b#]?', v_chord_root_str).group(0)

    print(f"This is a ii-V progression resolving to a new key.")
    print(f"The V chord is {v_chord_root_str}, and its root is {root_note_of_v}.")

    # 3. Find the index of this root note in our chromatic scale.
    try:
        # Handle cases where the note name might be a sharp, e.g., A# -> Bb
        if root_note_of_v == 'A#':
            root_note_of_v = 'Bb'
        elif root_note_of_v == 'C#':
            root_note_of_v = 'Db'
        elif root_note_of_v == 'D#':
            root_note_of_v = 'Eb'
        elif root_note_of_v == 'F#':
            root_note_of_v = 'Gb'
        elif root_note_of_v == 'G#':
            root_note_of_v = 'Ab'

        root_note_index = notes.index(root_note_of_v)
        print(f"The root note {root_note_of_v} is at index {root_note_index} in the chromatic scale.")
    except ValueError:
        print(f"Error: Could not find the root note {root_note_of_v} in the scale.")
        return

    # 4. A V chord resolves up a perfect fourth (5 semitones) to the tonic (I).
    # We calculate the index of the tonic note.
    semitones_in_perfect_fourth = 5
    tonic_index = (root_note_index + semitones_in_perfect_fourth) % 12
    
    print(f"The tonic note of the key is a perfect fourth (5 semitones) above the root of the V chord.")
    print(f"Calculation for the tonic's index: ({root_note_index} + {semitones_in_perfect_fourth}) % 12 = {tonic_index}")

    # 5. Find the tonic note using the new index.
    final_note = notes[tonic_index]
    print(f"The note at index {tonic_index} is {final_note}.")
    print("\nThe song's final phrase resolves to the key of " + final_note + ".")
    print("Therefore, the concluding note sung on the word 'you' is the tonic.")
    
    print(f"\nFinal Answer: {final_note}")

solve_final_note()