import bisect

def note_to_midi(note_str):
    """Converts a note string like 'C#4' to a MIDI number."""
    NOTES = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    # Handle cases like 'A#4' vs 'A4'
    if len(note_name) > 1 and note_name[1] in ('#', 'b'):
        pass
    else: # If note is natural, note_name is just the letter
        pass

    return (octave + 1) * 12 + NOTES[note_name]

def midi_to_note(midi_val):
    """Converts a MIDI number to a note string like 'C#4'."""
    NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    octave = (midi_val // 12) - 1
    note_index = midi_val % 12
    return f"{NOTE_NAMES[note_index]}{octave}"

def find_closest_dorian(note_val, dorian_notes_sorted):
    """Finds the closest note in the D Dorian scale."""
    # Find insertion point to locate neighbors
    idx = bisect.bisect_left(dorian_notes_sorted, note_val)

    # Exact match
    if idx < len(dorian_notes_sorted) and dorian_notes_sorted[idx] == note_val:
        return note_val

    # Handle edges
    if idx == 0:
        return dorian_notes_sorted[0]
    if idx == len(dorian_notes_sorted):
        return dorian_notes_sorted[-1]
    
    # Compare distances to neighbors
    low_note = dorian_notes_sorted[idx - 1]
    high_note = dorian_notes_sorted[idx]
    
    # If tie, choose the lower note
    if (note_val - low_note) <= (high_note - note_val):
        return low_note
    else:
        return high_note

def solve_music_sequence():
    """
    Applies a series of musical transformations to a sequence of notes.
    """
    # Initial sequence
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    sequence = [note_to_midi(n) for n in initial_notes_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    sequence = [n + 2 for n in sequence]

    # 2. Invert Around E4 (MIDI 64)
    axis_e4 = note_to_midi("E4")
    sequence = [(2 * axis_e4) - n for n in sequence]

    # 3. Retrograde the sequence
    sequence = sequence[::-1]

    # 4. Augment the intervals by adding 3 semitones to each interval
    sequence = [sequence[i] + 3 * i for i in range(len(sequence))]
    
    # 5. Change to Dorian Mode Starting from D4
    D_DORIAN_CLASSES = {0, 2, 4, 5, 7, 9, 11} # C, D, E, F, G, A, B
    d_dorian_midi_notes = sorted([i for i in range(128) if i % 12 in D_DORIAN_CLASSES])
    sequence = [find_closest_dorian(n, d_dorian_midi_notes) for n in sequence]

    # 6. Transpose down Minor Third (-3 semitones)
    sequence = [n - 3 for n in sequence]

    # 7. Invert around F4 (MIDI 65)
    axis_f4 = note_to_midi("F4")
    sequence = [(2 * axis_f4) - n for n in sequence]

    # 8. Transpose up one Octave (+12 semitones)
    sequence = [n + 12 for n in sequence]

    # Convert final MIDI sequence back to note names and print
    final_notes = [midi_to_note(n) for n in sequence]
    print(" ".join(final_notes))

solve_music_sequence()