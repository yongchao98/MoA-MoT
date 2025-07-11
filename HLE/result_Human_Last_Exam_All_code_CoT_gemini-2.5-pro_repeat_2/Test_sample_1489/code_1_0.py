import math

# --- Helper Functions for Note/MIDI Conversion ---

NOTE_TO_PC = {'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3, 'E': 4, 'F': 5, 'F#': 6, 'Gb': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11}
PC_TO_NOTE_SHARP = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

def note_to_midi(note_str):
    """Converts a note string like 'A4' to a MIDI number."""
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    pc = NOTE_TO_PC[note_name]
    return 12 * (octave + 1) + pc

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'A4'."""
    if not isinstance(midi_num, int):
        midi_num = int(round(midi_num))
    pc = midi_num % 12
    octave = midi_num // 12 - 1
    note_name = PC_TO_NOTE_SHARP[pc]
    return f"{note_name}{octave}"

def print_sequence(title, midi_sequence):
    """Prints the title and the sequence in note format."""
    note_sequence = [midi_to_note(m) for m in midi_sequence]
    print(f"{title:<45} {' '.join(note_sequence)}")

# --- Main Logic ---

# Initial Sequence
initial_notes = ['A4', 'C5', 'F4', 'E4', 'G4', 'C4', 'B4', 'D4']
midi_sequence = [note_to_midi(n) for n in initial_notes]
print_sequence("Initial Sequence:", midi_sequence)

# 1. Transpose up Major Second (+2 semitones)
midi_sequence = [n + 2 for n in midi_sequence]
print_sequence("1. Transpose up Major Second:", midi_sequence)

# 2. Invert Around E4 (MIDI 64)
axis = note_to_midi('E4')
midi_sequence = [2 * axis - n for n in midi_sequence]
print_sequence("2. Invert Around E4:", midi_sequence)

# 3. Retrograde the sequence
midi_sequence.reverse()
print_sequence("3. Retrograde:", midi_sequence)

# 4. Augment intervals by adding 3 semitones
intervals = [midi_sequence[i] - midi_sequence[i-1] for i in range(1, len(midi_sequence))]
augmented_intervals = [i + 3 for i in intervals]
new_sequence = [midi_sequence[0]]
for interval in augmented_intervals:
    new_sequence.append(new_sequence[-1] + interval)
midi_sequence = new_sequence
print_sequence("4. Augment Intervals by +3 semitones:", midi_sequence)

# 5. Change to Dorian Mode Starting from D4
d_dorian_pcs = {0, 2, 4, 5, 7, 9, 11}
quantized_sequence = []
for note in midi_sequence:
    # Find the nearest note in the D Dorian scale
    min_dist = float('inf')
    best_match = note
    # Search a range of octaves around the note
    for midi_val in range(note - 12, note + 12):
        if midi_val % 12 in d_dorian_pcs:
            dist = abs(note - midi_val)
            if dist < min_dist:
                min_dist = dist
                best_match = midi_val
            # Tie-breaking rule: if distances are equal, choose the lower note
            elif dist == min_dist:
                best_match = min(best_match, midi_val)
    quantized_sequence.append(best_match)
midi_sequence = quantized_sequence
print_sequence("5. Change to D Dorian Mode:", midi_sequence)

# 6. Transpose down Minor Third (-3 semitones)
midi_sequence = [n - 3 for n in midi_sequence]
print_sequence("6. Transpose down Minor Third:", midi_sequence)

# 7. Invert around F4 (MIDI 65)
axis = note_to_midi('F4')
midi_sequence = [2 * axis - n for n in midi_sequence]
print_sequence("7. Invert around F4:", midi_sequence)

# 8. Transposed up one Octave (+12 semitones)
midi_sequence = [n + 12 for n in midi_sequence]
print_sequence("8. Transpose up one Octave:", midi_sequence)

# --- Final Answer ---
final_notes = [midi_to_note(m) for m in midi_sequence]
final_answer_str = " ".join(final_notes)
print("\n" + "="*50)
print("Final sequence of notes:")
print(final_answer_str)
print("="*50)

# The final answer in the required format
# <<<A5 E6 C#5 F#5 B4 A4 C#5 G#4>>>