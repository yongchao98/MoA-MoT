import re

def note_to_midi(note_str):
    """Converts a note name (e.g., 'A#4') to a MIDI number."""
    NOTE_TO_MIDI_BASE = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    match = re.match(r'([A-G]#?)(-?\d+)', note_str.strip())
    if not match:
        raise ValueError(f"Invalid note format: {note_str}")
    name, octave = match.groups()
    octave = int(octave)
    return NOTE_TO_MIDI_BASE[name] + (octave + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number to a note name."""
    if not (0 <= midi_num <= 127):
        raise ValueError("MIDI number out of range (0-127)")
    MIDI_TO_NOTE_BASE = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}
    octave = (midi_num // 12) - 1
    note_index = midi_num % 12
    return f"{MIDI_TO_NOTE_BASE[note_index]}{octave}"

def transpose(midi_seq, semitones):
    """Transposes a sequence of MIDI notes by a number of semitones."""
    return [note + semitones for note in midi_seq]

def invert(midi_seq, axis_note_str):
    """Inverts a sequence of MIDI notes around a given axis note."""
    axis = note_to_midi(axis_note_str)
    return [2 * axis - note for note in midi_seq]

def retrograde(midi_seq):
    """Reverses a sequence of MIDI notes."""
    return midi_seq[::-1]

def augment_intervals(midi_seq, semitones):
    """Augments the intervals between consecutive notes."""
    if len(midi_seq) < 2:
        return midi_seq
    
    new_seq = [midi_seq[0]]
    for i in range(len(midi_seq) - 1):
        original_interval = midi_seq[i+1] - midi_seq[i]
        augmented_interval = original_interval + semitones
        new_note = new_seq[i] + augmented_interval
        new_seq.append(new_note)
    return new_seq

def change_to_mode(midi_seq, root_note_str, mode_pattern):
    """Maps notes to the closest notes in a given mode."""
    root_midi = note_to_midi(root_note_str)
    
    # Generate scale notes across a wide MIDI range
    scale_notes = set()
    for octave in range(-1, 10):
        for interval in mode_pattern:
            scale_notes.add(root_midi % 12 + interval + octave * 12)
            
    sorted_scale_notes = sorted(list(scale_notes))

    new_seq = []
    for note in midi_seq:
        # Find the minimum distance to a scale note
        min_dist = min(abs(s_note - note) for s_note in sorted_scale_notes)
        # Find all notes with that minimum distance
        closest_notes = [s_note for s_note in sorted_scale_notes if abs(s_note - note) == min_dist]
        # Choose the lowest note in case of a tie
        new_seq.append(min(closest_notes))
        
    return new_seq

# Initial sequence
initial_notes = "A4 C5 F4 E4 G4 C4 B4 D4".split()

# Convert to MIDI
current_sequence = [note_to_midi(n) for n in initial_notes]

# 1. Transpose up Major Second (+2)
current_sequence = transpose(current_sequence, 2)

# 2. Invert Around E4
current_sequence = invert(current_sequence, 'E4')

# 3. Retrograde the sequence
current_sequence = retrograde(current_sequence)

# 4. Augment the intervals by adding 3 semitones
current_sequence = augment_intervals(current_sequence, 3)

# 5. Change to Dorian Mode Starting from D4
dorian_pattern = [0, 2, 3, 5, 7, 9, 10] # T, T, S, T, T, T, S intervals from root -> semitones 0, 2, 3, 5, 7, 9, 10
current_sequence = change_to_mode(current_sequence, 'D4', dorian_pattern)

# 6. Transpose down Minor Third (-3)
current_sequence = transpose(current_sequence, -3)

# 7. Invert around F4
current_sequence = invert(current_sequence, 'F4')

# 8. Transposed up one Octave (+12)
current_sequence = transpose(current_sequence, 12)

# Convert final MIDI sequence back to note names
final_notes = [midi_to_note(n) for n in current_sequence]

# Print the final result
print(" ".join(final_notes))