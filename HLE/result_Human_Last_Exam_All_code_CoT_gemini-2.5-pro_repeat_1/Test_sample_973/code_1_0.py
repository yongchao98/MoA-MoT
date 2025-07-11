def solve_music_puzzle():
    """
    This script solves the music theory puzzle based on the provided chord progression.
    """
    
    # The chromatic scale using flats, which are common in jazz progressions.
    # Index: 0   1    2   3    4   5   6    7   8    9   10   11
    notes = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B']
    
    # The final chord pair for "bir-day" is Cm7 - F7(9).
    # This is a ii-V progression. We focus on the V chord, F7(9).
    v_chord_root = 'F'
    
    print("Step 1: Identify the final harmonic progression for '...happy birthday...'.")
    print("The chords are Cm7 - F7(9). This is a 'ii-V' progression.")
    print("-" * 20)
    
    print("Step 2: Determine the implied key from the 'V' chord.")
    print("The 'V' chord is F7(9), which has a root note of F.")
    
    # In music theory, the tonic (I) is a perfect fifth (7 semitones) below the root of the V chord.
    semitones_in_perfect_fifth = 7
    
    # Find the index of the V chord's root.
    v_root_index = notes.index(v_chord_root)
    
    # Calculate the index of the tonic note.
    # We subtract 7 semitones and use modulo to wrap around the scale.
    # (v_root_index - 7 + 12) % 12 ensures the result is positive.
    tonic_index = (v_root_index - semitones_in_perfect_fifth + len(notes)) % len(notes)
    
    # Get the tonic note from the list.
    tonic_note = notes[tonic_index]
    
    print(f"The tonic of the key is a perfect fifth ({semitones_in_perfect_fifth} semitones) below {v_chord_root}.")
    print(f"Calculation: The note {semitones_in_perfect_fifth} semitones below {v_chord_root} is {tonic_note}.")
    print("-" * 20)
    
    print("Step 3: Identify the final melodic note.")
    print("The melody of 'Happy Birthday to You' ends on the tonic of the key.")
    print("The final note, sung on the word 'you', is therefore the tonic.")
    print("-" * 20)
    
    print(f"Final Answer: The note is {tonic_note}.")

solve_music_puzzle()
<<<Bb>>>