import math

# --- Helper Functions ---

NOTE_NAMES_SHARP = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
NOTE_MAP = {name: i for i, name in enumerate(NOTE_NAMES_SHARP)}

def note_to_midi(note_str):
    """Converts a note name (e.g., 'A4', 'C#5') to a MIDI number."""
    note_str = note_str.upper()
    octave = int(note_str[-1])
    name = note_str[:-1]
    if name not in NOTE_MAP:
        raise ValueError(f"Unknown note name: {name}")
    note_index = NOTE_MAP[name]
    return 12 * (octave + 1) + note_index

def midi_to_note(midi_num):
    """Converts a MIDI number to a note name (e.g., 'C#5')."""
    if not isinstance(midi_num, int):
        raise TypeError("MIDI number must be an integer.")
    octave = midi_num // 12 - 1
    note_index = midi_num % 12
    return f"{NOTE_NAMES_SHARP[note_index]}{octave}"

# --- Transformation Functions ---

def transpose(midi_notes, semitones):
    """Transposes a sequence of MIDI notes by a number of semitones."""
    return [note + semitones for note in midi_notes]

def invert(midi_notes, axis_note_str):
    """Inverts a sequence of MIDI notes around an axis note."""
    axis_midi = note_to_midi(axis_note_str)
    return [2 * axis_midi - note for note in midi_notes]

def retrograde(midi_notes):
    """Reverses the order of a sequence of MIDI notes."""
    return midi_notes[::-1]

def augment_intervals(midi_notes, semitone_addition):
    """Augments the intervals between consecutive notes."""
    if len(midi_notes) < 2:
        return midi_notes
    new_notes = [midi_notes[0]]
    for i in range(len(midi_notes) - 1):
        interval = midi_notes[i+1] - midi_notes[i]
        augmented_interval = interval + semitone_addition
        next_note = new_notes[i] + augmented_interval
        new_notes.append(next_note)
    return new_notes

def snap_to_mode(midi_notes, root_note_str, mode_intervals):
    """Snaps notes to the closest note in a specified mode."""
    root_midi = note_to_midi(root_note_str)
    scale_notes = set()
    # Generate scale notes across several octaves
    for oct_shift in range(-2, 4):
        current_note = root_midi + (12 * oct_shift)
        scale_notes.add(current_note)
        for interval in mode_intervals:
            current_note += interval
            scale_notes.add(current_note)

    sorted_scale = sorted(list(scale_notes))
    snapped_notes = []
    for note in midi_notes:
        # Find the closest note in the scale.
        # In case of a tie, choose the higher note.
        min_diff = float('inf')
        closest_notes = []
        for scale_note in sorted_scale:
            diff = abs(note - scale_note)
            if diff < min_diff:
                min_diff = diff
                closest_notes = [scale_note]
            elif diff == min_diff:
                closest_notes.append(scale_note)
        snapped_notes.append(max(closest_notes))
    return snapped_notes

# --- Main Execution ---

# Initial sequence
initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
sequence = [note_to_midi(n) for n in initial_notes_str.split()]

# 1. Transpose up Major Second (+2 semitones)
sequence = transpose(sequence, 2)

# 2. Invert Around E4
sequence = invert(sequence, "E4")

# 3. Retrograde the sequence
sequence = retrograde(sequence)

# 4. Augment the intervals by adding 3 semitones
sequence = augment_intervals(sequence, 3)

# 5. Change to Dorian Mode Starting from D4
# Dorian mode intervals: Tone, Semi-tone, Tone, Tone, Tone, Semi-tone, Tone
dorian_intervals = [2, 1, 2, 2, 2, 1, 2]
sequence = snap_to_mode(sequence, "D4", dorian_intervals)

# 6. Transpose down Minor Third (-3 semitones)
sequence = transpose(sequence, -3)

# 7. Invert around F4
sequence = invert(sequence, "F4")

# 8. Transposed up one Octave (+12 semitones)
sequence = transpose(sequence, 12)

# Convert final MIDI sequence back to note names and print
final_note_names = [midi_to_note(n) for n in sequence]
print(' '.join(final_note_names))
