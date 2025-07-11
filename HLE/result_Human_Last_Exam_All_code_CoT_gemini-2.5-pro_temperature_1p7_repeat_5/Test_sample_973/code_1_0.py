def solve_birthday_note():
    """
    This script solves the musical puzzle by determining the song's key
    and applying a fundamental rule of the "Happy Birthday" melody.
    """

    # --- Step 1: Analyze the chord progression pattern ---

    # The progression for the first "birthday" and "to you"
    first_birthday_chords = ("Cm7", "F7(9)")
    first_you_chord = "Bm7"

    # The progression for the final "birthday"
    final_birthday_chords = ("Cm7", "F7(9)")

    # The problem states the pattern is consistent.
    # We deduce the chord for the final "you" by applying the pattern.
    final_you_chord = ""
    if final_birthday_chords == first_birthday_chords:
        final_you_chord = first_you_chord

    print(f"The chords for the final 'birthday' are: {final_birthday_chords[0]} {final_birthday_chords[1]}")
    print(f"Based on the song's established pattern, the chord for the final 'you' is: {final_you_chord}")

    # --- Step 2: Determine the key from the primary ii-V progression ---

    ii_chord_root = final_birthday_chords[0].replace('m7', '') # Gets 'C' from 'Cm7'

    # The key (I chord) of a ii-V-I progression is a whole step (2 semitones) below the ii chord's root.
    # We'll use a chromatic scale to calculate this.
    notes_sharp = ['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#']
    
    ii_root_index = notes_sharp.index(ii_chord_root)
    tonic_index = (ii_root_index - 2 + 12) % 12
    
    # The note at the tonic_index is A#, which is commonly written as Bb.
    song_key = "Bb"

    print(f"The progression {final_birthday_chords[0]}-{final_birthday_chords[1]} is a 'ii-V', which establishes the song's key as: {song_key}")

    # --- Step 3: Apply the melodic rule ---
    
    # The traditional melody of "Happy Birthday" concludes on the tonic of the key.
    final_melody_note = song_key
    
    print("The melody of 'Happy Birthday' traditionally ends on the tonic (root note) of the key.")
    print(f"Therefore, the note sung on the concluding word 'you' is: {final_melody_note}")


solve_birthday_note()
<<<Bb>>>