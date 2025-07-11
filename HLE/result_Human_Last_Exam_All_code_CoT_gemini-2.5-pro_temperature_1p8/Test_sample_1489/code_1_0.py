import math

def notes_to_midi(note_strs):
    """Converts a space-separated string of notes (e.g., 'A4 C#5') to a list of MIDI numbers."""
    NOTE_VALS = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    midi_nums = []
    for note_str in note_strs.split():
        name = note_str[0].upper()
        octave_part = note_str[1:]
        
        accidental = 0
        if '#' in octave_part:
            accidental = 1
            octave_part = octave_part.replace('#', '')
        elif 'b' in octave_part:
            accidental = -1
            octave_part = octave_part.replace('b', '')
            
        octave = int(octave_part)
        note_val = NOTE_VALS[name]
        
        midi_num = 12 * (octave + 1) + note_val + accidental
        midi_nums.append(midi_num)
    return midi_nums

def midi_to_notes(midi_nums):
    """Converts a list of MIDI numbers to a space-separated string of note names."""
    NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    note_strs = []
    for midi_num in midi_nums:
        note_name = NOTE_NAMES[midi_num % 12]
        octave = (midi_num // 12) - 1
        note_strs.append(f"{note_name}{octave}")
    return " ".join(note_strs)

# 1. Initial sequence
initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
notes = notes_to_midi(initial_notes_str)

# Operation 1: Transpose up Major Second (+2 semitones)
notes = [n + 2 for n in notes]

# Operation 2: Invert Around E4 (MIDI 64)
axis_e4 = 64
notes = [(2 * axis_e4) - n for n in notes]

# Operation 3: Retrograde the sequence
notes = notes[::-1]

# Operation 4: Augment the intervals by adding 3 semitones
if len(notes) > 1:
    intervals = [notes[i+1] - notes[i] for i in range(len(notes)-1)]
    new_intervals = [i + 3 for i in intervals]
    new_notes = [notes[0]]
    for interval in new_intervals:
        new_notes.append(new_notes[-1] + interval)
    notes = new_notes

# Operation 5: Change to Dorian Mode Starting from D4
root_d4_midi = notes_to_midi("D4")[0]
dorian_intervals = [2, 1, 2, 2, 2, 1, 2] # T-S-T-T-T-S-T
dorian_degrees = {0}
current_degree = 0
for interval in dorian_intervals:
    current_degree += interval
    dorian_degrees.add(current_degree % 12)

# Generate a wide range of scale notes
scale_notes = []
for octave in range(-1, 10):
    for degree in dorian_degrees:
        scale_notes.append(root_d4_midi % 12 + 12 * octave + degree)
scale_notes.sort()

# Quantize notes to the scale
quantized_notes = []
for note in notes:
    min_dist = float('inf')
    candidates = []
    for scale_note in scale_notes:
        dist = abs(note - scale_note)
        if dist < min_dist:
            min_dist = dist
            candidates = [scale_note]
        elif dist == min_dist:
            candidates.append(scale_note)
    # Tie-break by choosing the lower note
    quantized_notes.append(min(candidates))

# Transpose so the sequence starts on D4
tonic_d4_midi = notes_to_midi("D4")[0]
if quantized_notes:
    transposition = tonic_d4_midi - quantized_notes[0]
    notes = [n + transposition for n in quantized_notes]

# Operation 6: Transpose down Minor Third (-3 semitones)
notes = [n - 3 for n in notes]

# Operation 7: Invert around F4 (MIDI 65)
axis_f4 = 65
notes = [(2 * axis_f4) - n for n in notes]

# Operation 8: Transposed up one Octave (+12 semitones)
notes = [n + 12 for n in notes]

# Final Answer
final_note_sequence = midi_to_notes(notes)
print(final_note_sequence)