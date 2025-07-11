import math

def note_to_midi(note_str):
    """Converts a note string like 'A#4' or 'F4' to its MIDI number."""
    NOTE_MAP = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    # Handle two-character note names like 'C#'
    if len(note_name) > 1 and note_name[1] in ['#', 'b']:
        pass
    else: # Handle single-character note names like 'A', 'C'
        note_name = note_str[0]
        octave = int(note_str[1:])

    # The formula MIDI = 12 * (octave + 1) + note_index is a standard
    return NOTE_MAP[note_name] + (octave + 1) * 12

def midi_to_note(midi_val):
    """Converts a MIDI number to its note string representation, using sharps."""
    NOTE_MAP_SHARP = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}
    octave = (midi_val // 12) - 1
    note_name = NOTE_MAP_SHARP[midi_val % 12]
    return f"{note_name}{octave}"

def is_in_scale(note, root, scale_intervals):
    """Checks if a MIDI note is part of a given scale."""
    return (note - root) % 12 in scale_intervals

def find_closest_in_scale(note, root, scale_intervals):
    """Finds the closest note in a scale, rounding up on ties."""
    if is_in_scale(note, root, scale_intervals):
        return note
    
    delta = 1
    while True:
        lower = note - delta
        upper = note + delta
        
        lower_in_scale = is_in_scale(lower, root, scale_intervals)
        upper_in_scale = is_in_scale(upper, root, scale_intervals)

        # On a tie, the convention is to round up to the higher pitch.
        if lower_in_scale and upper_in_scale:
            return upper
        elif lower_in_scale:
            return lower
        elif upper_in_scale:
            return upper
        
        delta += 1

def solve_music_sequence():
    """
    Solves the musical transformation problem step-by-step.
    """
    # Initial sequence
    initial_note_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    
    # Step 0: Convert initial sequence to MIDI numbers
    midi_seq = [note_to_midi(n) for n in initial_note_str.split()]
    
    # Operation 1: Transpose up Major Second (+2 semitones)
    midi_seq = [n + 2 for n in midi_seq]

    # Operation 2: Invert Around E4
    pivot_e4 = note_to_midi('E4')
    midi_seq = [2 * pivot_e4 - n for n in midi_seq]
    
    # Operation 3: Retrograde the sequence
    midi_seq.reverse()

    # Operation 4: Augment intervals by adding 3 semitones
    if len(midi_seq) > 1:
        augmented_seq = [midi_seq[0]]
        for i in range(len(midi_seq) - 1):
            interval = midi_seq[i+1] - midi_seq[i]
            new_interval = interval + 3
            augmented_seq.append(augmented_seq[i] + new_interval)
        midi_seq = augmented_seq

    # Operation 5: Change to Dorian Mode Starting from D4
    root_d4 = note_to_midi('D4')
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10} # W-H-W-W-W-H-W
    midi_seq = [find_closest_in_scale(n, root_d4, dorian_intervals) for n in midi_seq]

    # Operation 6: Transpose down Minor Third (-3 semitones)
    midi_seq = [n - 3 for n in midi_seq]

    # Operation 7: Invert around F4
    pivot_f4 = note_to_midi('F4')
    midi_seq = [2 * pivot_f4 - n for n in midi_seq]

    # Operation 8: Transposed up one Octave (+12 semitones)
    midi_seq = [n + 12 for n in midi_seq]

    # Final step: Convert back to note names and print
    final_notes = [midi_to_note(n) for n in midi_seq]
    print(" ".join(final_notes))

# Execute the solution
solve_music_sequence()