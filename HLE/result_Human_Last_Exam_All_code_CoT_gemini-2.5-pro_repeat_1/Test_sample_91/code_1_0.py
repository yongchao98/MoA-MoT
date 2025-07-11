def find_matching_scales():
    """
    Identifies possible musical scales that contain a given set of notes.
    The user should replace the example notes with their own transcription.
    """
    # --- Step 1: Define notes and scales ---

    # Helper dictionaries to map between note names and numerical values (0-11)
    NOTE_TO_VAL = {
        'C': 0, 'C#': 1, 'DB': 1, 'D': 2, 'D#': 3, 'EB': 3, 'E': 4,
        'F': 5, 'F#': 6, 'GB': 6, 'G': 7, 'G#': 8, 'AB': 8, 'A': 9,
        'A#': 10, 'BB': 10, 'B': 11
    }
    VAL_TO_NOTE = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    # Define common scales by their intervals in semitones from the root note
    SCALES = {
        'Major': [0, 2, 4, 5, 7, 9, 11],
        'Natural Minor': [0, 2, 3, 5, 7, 8, 10],
        'Harmonic Minor': [0, 2, 3, 5, 7, 8, 11],
        'Major Pentatonic': [0, 2, 4, 7, 9],
        'Minor Pentatonic': [0, 3, 5, 7, 10],
        'Blues': [0, 3, 5, 6, 7, 10],
        'Whole Tone': [0, 2, 4, 6, 8, 10],
        'Diminished (WH)': [0, 2, 3, 5, 6, 8, 9, 11]
    }

    # --- Step 2: Input the transcribed notes ---

    # !!! IMPORTANT !!!
    # Replace the notes in this list with the melody notes you hear in the song
    # between 0:39 and 0:43. This is just a placeholder example.
    transcribed_notes = ['G', 'A', 'B', 'C#', 'D#']

    print(f"Analyzing notes: {transcribed_notes}")

    # --- Step 3: Analyze and find matching scales ---

    try:
        # Convert the input notes to a set of unique numerical values
        note_values = {NOTE_TO_VAL[note.upper()] for note in transcribed_notes}
    except KeyError as e:
        print(f"Error: Unrecognized note name {e}. Please use standard notation (e.g., 'C', 'F#', 'Bb').")
        return

    # Print the numerical representation of the notes, as requested.
    # An "equation" isn't applicable here, so we print the component numbers.
    print(f"Numerical representation (where C=0): {sorted(list(note_values))}\n")

    found_matches = []
    # Iterate through all 12 possible root notes (C, C#, D, etc.)
    for root_val in range(12):
        # Iterate through all the scale types we defined
        for scale_name, intervals in SCALES.items():
            # Construct the full set of notes for the current scale and root
            current_scale_notes = {(root_val + interval) % 12 for interval in intervals}
            
            # Check if all the transcribed notes are present in the current scale
            if note_values.issubset(current_scale_notes):
                root_note_name = VAL_TO_NOTE[root_val]
                found_matches.append(f"{root_note_name} {scale_name}")

    # --- Step 4: Print the results ---
    if found_matches:
        print("The provided notes fit into the following possible scales:")
        for match in found_matches:
            print(f"- {match}")
    else:
        print("No common scales were found containing all the provided notes.")
        print("The passage might be chromatic or use a more exotic scale not listed.")

# Run the analysis
find_matching_scales()
