import math

# --- Helper Functions for Note/MIDI Conversion ---

NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
NOTE_MAP = {name: i for i, name in enumerate(NOTE_NAMES)}

def note_to_midi(note_str):
    """Converts a note string like 'A#4' to a MIDI number."""
    note_str = note_str.upper()
    if len(note_str) > 2 and (note_str[1] == '#' or note_str[1] == 'B'):
        name = note_str[:2]
        octave = int(note_str[2:])
    else:
        name = note_str[:1]
        octave = int(note_str[1:])
    
    # Using 'B' for flat is not standard here, so we assume '#' notation
    if 'B' in name and name != 'B':
        # Handle flats by converting to sharps, e.g., Db -> C#
        base_note_index = NOTE_MAP[name[0]]
        pitch_class = (base_note_index - 1 + 12) % 12
    else:
        pitch_class = NOTE_MAP[name]
        
    return (octave + 1) * 12 + pitch_class

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'A#4'."""
    if not isinstance(midi_num, int):
        return "Invalid"
    octave = (midi_num // 12) - 1
    pitch_class = midi_num % 12
    return f"{NOTE_NAMES[pitch_class]}{octave}"

# --- Main Logic ---

def solve_music_sequence():
    """
    Applies a series of musical operations to a sequence of notes.
    """
    # Initial sequence
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    notes = [note_to_midi(n) for n in initial_notes_str]
    # Initial MIDI: [69, 72, 65, 64, 67, 60, 71, 62]

    # 1. Transpose up Major Second (2 semitones)
    notes = [n + 2 for n in notes]
    # Result: [71, 74, 67, 66, 69, 62, 73, 64]

    # 2. Invert Around E4 (MIDI 64)
    axis_e4 = note_to_midi('E4')
    notes = [(2 * axis_e4) - n for n in notes]
    # Result: [57, 54, 61, 62, 59, 66, 55, 64]

    # 3. Retrograde the sequence
    notes = notes[::-1]
    # Result: [64, 55, 66, 59, 62, 61, 54, 57]

    # 4. Augment intervals by adding 3 semitones
    original_intervals = [notes[i] - notes[i-1] for i in range(1, len(notes))]
    augmented_intervals = [i + 3 for i in original_intervals]
    augmented_notes = [notes[0]]
    for interval in augmented_intervals:
        augmented_notes.append(augmented_notes[-1] + interval)
    notes = augmented_notes
    # Result: [64, 58, 72, 68, 74, 76, 72, 78]

    # 5. Change to Dorian Mode Starting from D4
    dorian_root = note_to_midi('D4')
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10} # Semitones from root for Dorian
    dorian_pitch_classes = {(dorian_root % 12 + i) % 12 for i in dorian_intervals}
    
    mapped_notes = []
    for n in notes:
        pitch_class = n % 12
        if pitch_class in dorian_pitch_classes:
            mapped_notes.append(n)
        else:
            # Find the nearest scale tone, round up on ties
            d_down = 1
            while (pitch_class - d_down + 12) % 12 not in dorian_pitch_classes:
                d_down += 1
            
            d_up = 1
            while (pitch_class + d_up) % 12 not in dorian_pitch_classes:
                d_up += 1

            if d_up < d_down:
                mapped_notes.append(n + d_up)
            elif d_down < d_up:
                mapped_notes.append(n - d_down)
            else: # Tie, round up
                mapped_notes.append(n + d_up)
    notes = mapped_notes
    # Result: [64, 59, 72, 69, 74, 76, 72, 79]

    # 6. Transpose down Minor Third (3 semitones)
    notes = [n - 3 for n in notes]
    # Result: [61, 56, 69, 66, 71, 73, 69, 76]

    # 7. Invert around F4 (MIDI 65)
    axis_f4 = note_to_midi('F4')
    notes = [(2 * axis_f4) - n for n in notes]
    # Result: [69, 74, 61, 64, 59, 57, 61, 54]

    # 8. Transposed up one Octave (12 semitones)
    notes = [n + 12 for n in notes]
    # Result: [81, 86, 73, 76, 71, 69, 73, 66]

    # Convert final MIDI numbers back to note strings
    final_note_sequence = [midi_to_note(n) for n in notes]
    
    # Print the final result
    print(" ".join(final_note_sequence))

solve_music_sequence()