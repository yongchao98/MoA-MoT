import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
    notes_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    note = note_str[0].upper()
    octave_part = note_str[1:]
    
    is_sharp = '#' in octave_part
    octave = int(octave_part.replace('#', ''))
    
    midi_val = 12 * (octave + 1) + notes_map[note]
    if is_sharp:
        midi_val += 1
        
    return midi_val

def midi_to_note(midi_num):
    """Converts a MIDI number back to its note string representation."""
    semitone_map = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    octave = (midi_num // 12) - 1
    note_name = semitone_map[midi_num % 12]
    return f"{note_name}{octave}"

def quantize_to_d_dorian(midi_note):
    """Finds the closest note in the D Dorian scale to a given MIDI note."""
    # D Dorian scale tones (pitch classes, where C=0): D(2), E(4), F(5), G(7), A(9), B(11), C(0)
    scale_degrees = {0, 2, 4, 5, 7, 9, 11}
    
    # If the note is already in the scale, return it
    if midi_note % 12 in scale_degrees:
        return midi_note

    # Find the nearest note above and below that are in the scale
    dist_down = 1
    while (midi_note - dist_down) % 12 not in scale_degrees:
        dist_down += 1
    
    dist_up = 1
    while (midi_note + dist_up) % 12 not in scale_degrees:
        dist_up += 1

    # Return the closest note, favoring the lower note in case of a tie
    if dist_down <= dist_up:
        return midi_note - dist_down
    else:
        return midi_note + dist_up

def solve_music_sequence():
    """
    Performs the eight specified musical operations on the initial sequence of notes.
    """
    # Initial sequence
    initial_notes = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # Operation 1: Transpose up Major Second (+2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]

    # Operation 2: Invert Around E4
    axis_e4 = note_to_midi('E4')
    midi_sequence = [(2 * axis_e4) - n for n in midi_sequence]

    # Operation 3: Retrograde the sequence
    midi_sequence.reverse()

    # Operation 4: Augment the intervals by adding 3 semitones
    augmented_sequence = [midi_sequence[0]]
    for i in range(1, len(midi_sequence)):
        interval = midi_sequence[i] - midi_sequence[i-1]
        new_interval = interval + 3
        next_note = augmented_sequence[-1] + new_interval
        augmented_sequence.append(next_note)
    midi_sequence = augmented_sequence

    # Operation 5: Change to Dorian Mode Starting from D4
    midi_sequence = [quantize_to_d_dorian(n) for n in midi_sequence]

    # Operation 6: Transpose down Minor Third (-3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]
    
    # Operation 7: Invert around F4
    axis_f4 = note_to_midi('F4')
    midi_sequence = [(2 * axis_f4) - n for n in midi_sequence]
    
    # Operation 8: Transpose up one Octave (+12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]

    # Convert the final MIDI numbers back to note names
    final_notes = [midi_to_note(n) for n in midi_sequence]

    # Print the final result
    print(" ".join(final_notes))

solve_music_sequence()
<<<A5 E6 C#5 F#5 B4 A4 C#5 G#4>>>