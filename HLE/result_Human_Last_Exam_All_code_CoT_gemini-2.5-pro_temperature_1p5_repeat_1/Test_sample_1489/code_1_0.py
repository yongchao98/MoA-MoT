import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C5#' into a MIDI number."""
    note_map = {'C': 0, 'C#': 1, 'DB': 1, 'D': 2, 'D#': 3, 'EB': 3, 'E': 4, 'F': 5, 
                'F#': 6, 'GB': 6, 'G': 7, 'G#': 8, 'AB': 8, 'A': 9, 'A#': 10, 'BB': 10, 'B': 11}
    
    note_name = note_str[:-1]
    octave_char = note_str[-1]

    if note_name[-1] in ['#', 'B']:
        accidental = note_name[-1]
        base_note = note_name[:-1]
        if accidental == 'B':  # handle flats
             note_name = base_note + 'b'.upper()
        else:
             note_name = base_note + '#'
        octave_char = note_str[-1]
        octave = int(octave_char)
    else:
        note_name = note_str[:-1]
        octave = int(note_str[-1])
        
    # Handle cases like C#5 -> C#, 5
    if len(note_name) > 1 and note_name[-1] in ['#', 'b']:
        base_note_name = note_name
    else: # e.g. A4
        base_note_name = note_name

    semitone = note_map[base_note_name.upper()]
    return (octave + 1) * 12 + semitone

def midi_to_note(midi):
    """Converts a MIDI number into a note string like 'A4' or 'C5#'."""
    note_names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    octave = midi // 12 - 1
    note_index = midi % 12
    base_name = note_names[note_index]
    
    if len(base_name) > 1: # It's a sharp note
        return f"{base_name[0]}{octave}{base_name[1]}"
    else: # It's a natural note
        return f"{base_name}{octave}"

def main():
    # --- Initial Sequence ---
    initial_note_strs = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    notes = [note_to_midi(n) for n in initial_note_strs]
    print(f"Initial MIDI sequence:                     {notes}")

    # --- Operation 1: Transpose up Major Second (+2 semitones) ---
    notes = [n + 2 for n in notes]
    print(f"1. After Transposing up a Major Second:    {notes}")

    # --- Operation 2: Invert Around E4 ---
    inversion_axis_e4 = note_to_midi("E4")
    notes = [(2 * inversion_axis_e4) - n for n in notes]
    print(f"2. After Inverting around E4 (MIDI {inversion_axis_e4}):      {notes}")

    # --- Operation 3: Retrograde ---
    notes.reverse()
    print(f"3. After Retrograde:                       {notes}")

    # --- Operation 4: Augment intervals by adding 3 semitones ---
    intervals = [notes[i] - notes[i - 1] for i in range(1, len(notes))]
    augmented_intervals = [i + 3 for i in intervals]
    augmented_notes = [notes[0]]
    for interval in augmented_intervals:
        augmented_notes.append(augmented_notes[-1] + interval)
    notes = augmented_notes
    print(f"4. After Augmenting Intervals by 3:        {notes}")

    # --- Operation 5: Change to Dorian Mode Starting from D4 ---
    d_dorian_root = note_to_midi("D4")
    d_dorian_intervals = [0, 2, 3, 5, 7, 9, 10]  # W-H-W-W-W-H-W intervals in semitones
    d_dorian_scale_notes = set()
    for o in range(-2, 10): # Populate scale across several octaves
        for i in d_dorian_intervals:
            d_dorian_scale_notes.add(d_dorian_root + o * 12 + i - (d_dorian_root%12))

    # Transpose sequence so the first note is D4
    transposition_val = d_dorian_root - notes[0]
    transposed_notes = [n + transposition_val for n in notes]
    
    # Quantize to the nearest Dorian note
    quantized_notes = []
    for note in transposed_notes:
        # Find the nearest Dorian note. Ties are broken by rounding up.
        diffs = {abs(note - dorian_note): dorian_note for dorian_note in d_dorian_scale_notes}
        min_diff = min(diffs.keys())
        
        # Check for ties
        candidates = [v for k, v in diffs.items() if k == min_diff]
        quantized_notes.append(max(candidates)) # Round up on a tie
        
    notes = quantized_notes
    print(f"5. After Changing to D Dorian from D4:     {notes}")

    # --- Operation 6: Transpose down Minor Third (-3 semitones) ---
    notes = [n - 3 for n in notes]
    print(f"6. After Transposing down a Minor Third:   {notes}")

    # --- Operation 7: Invert around F4 ---
    inversion_axis_f4 = note_to_midi("F4")
    notes = [(2 * inversion_axis_f4) - n for n in notes]
    print(f"7. After Inverting around F4 (MIDI {inversion_axis_f4}):      {notes}")
    
    # --- Operation 8: Transpose up one Octave (+12 semitones) ---
    notes = [n + 12 for n in notes]
    print(f"8. After Transposing up an Octave:         {notes}")

    # --- Final Conversion and Output ---
    final_note_sequence = " ".join([midi_to_note(n) for n in notes])
    print("\nFinal note sequence:")
    print(final_note_sequence)

if __name__ == "__main__":
    main()

<<<B5 E6 D5 G5# C5# B4 D5 A4>>>