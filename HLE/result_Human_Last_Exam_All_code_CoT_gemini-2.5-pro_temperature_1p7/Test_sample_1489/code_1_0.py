import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C#5' to a MIDI number."""
    NOTES = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    note_name = note_str[0].upper()
    octave_str = note_str[1:]
    
    accidental = 0
    if len(octave_str) > 1 and octave_str[0] in '#b':
        if octave_str[0] == '#':
            accidental = 1
        elif octave_str[0] == 'b':
            accidental = -1
        octave_str = octave_str[1:]
        
    octave = int(octave_str)
    pitch_class = NOTES[note_name] + accidental
    return 12 * (octave + 1) + pitch_class

def midi_to_note(midi_val):
    """Converts a MIDI number to a note string like 'A4' or 'C#5'."""
    NOTES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    octave = (midi_val // 12) - 1
    note_index = midi_val % 12
    note_name = NOTES[note_index]
    return f"{note_name}{octave}"

def apply_operations(initial_sequence):
    """Applies the eight musical operations to the initial sequence."""
    
    # Initial sequence
    midi_seq = [note_to_midi(n) for n in initial_sequence.split()]
    # [69, 72, 65, 64, 67, 60, 71, 62]

    # 1. Transpose up Major Second (2 semitones)
    midi_seq = [n + 2 for n in midi_seq]
    # [71, 74, 67, 66, 69, 62, 73, 64]

    # 2. Invert Around E4
    axis_e4 = note_to_midi('E4') # 64
    midi_seq = [(2 * axis_e4) - n for n in midi_seq]
    # [57, 54, 61, 62, 59, 66, 55, 64]

    # 3. Retrograde the sequence
    midi_seq = midi_seq[::-1]
    # [64, 55, 66, 59, 62, 61, 54, 57]
    
    # 4. Augment intervals by adding 3 semitones
    new_seq = [midi_seq[0]]
    for i in range(len(midi_seq) - 1):
        interval = midi_seq[i+1] - midi_seq[i]
        new_interval = interval + 3
        next_note = new_seq[-1] + new_interval
        new_seq.append(next_note)
    midi_seq = new_seq
    # [64, 58, 72, 68, 74, 76, 72, 78]

    # 5. Change to Dorian Mode Starting from D4
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10}
    root_d4 = note_to_midi('D4') # 62
    
    # Generate a list of all notes in the D Dorian scale over a wide range
    dorian_scale_notes = []
    for midi_note in range(128):
        if (midi_note - root_d4) % 12 in dorian_intervals:
            dorian_scale_notes.append(midi_note)

    new_seq = []
    for note in midi_seq:
        # Find the closest note in the scale. In case of a tie, min returns the first one (lower pitch).
        closest_note = min(dorian_scale_notes, key=lambda x: abs(x - note))
        new_seq.append(closest_note)
    midi_seq = new_seq
    # [64, 57, 72, 67, 74, 76, 72, 77]
    
    # 6. Transpose down Minor Third (3 semitones)
    midi_seq = [n - 3 for n in midi_seq]
    # [61, 54, 69, 64, 71, 73, 69, 74]

    # 7. Invert around F4
    axis_f4 = note_to_midi('F4') # 65
    midi_seq = [(2 * axis_f4) - n for n in midi_seq]
    # [69, 76, 61, 66, 59, 57, 61, 56]

    # 8. Transposed up one Octave (12 semitones)
    midi_seq = [n + 12 for n in midi_seq]
    # [81, 88, 73, 78, 71, 69, 73, 68]
    
    # Convert final MIDI sequence back to note names
    final_notes = [midi_to_note(n) for n in midi_seq]
    return " ".join(final_notes)

# The initial sequence of eight musical notes
initial_notes_sequence = "A4 C5 F4 E4 G4 C4 B4 D4"

# Calculate the final sequence
final_sequence = apply_operations(initial_notes_sequence)

# Print the final answer
print(final_sequence)
<<<A5 E6 C5# F5# B4 A4 C5# G4#>>>