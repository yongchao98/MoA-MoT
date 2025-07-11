import math

def note_to_midi(note_str):
    """Converts a note string like 'A4' or 'C5#' to a MIDI number."""
    NOTES = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    note_name = note_str[:-1]
    octave = int(note_str[-1])
    
    if len(note_name) > 1 and (note_name[1] == '#' or note_name[1] == 'b'):
        # Handles cases like 'C#4'
        note = note_name[:2]
        octave = int(note_str[2:])
    else:
        # Handles cases like 'A4'
        note = note_name[0]
        octave = int(note_str[1:])
        
    semitone = NOTES[note]
    # MIDI formula: C4 = 60
    return 12 * (octave + 1) + semitone

def midi_to_note(midi_num):
    """Converts a MIDI number to a note string like 'A4' or 'C#5'."""
    NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    
    octave = midi_num // 12 - 1
    note_index = midi_num % 12
    note_name = NOTE_NAMES[note_index]
    
    return f"{note_name}{octave}"

def get_dorian_scale_notes(start_midi, num_octaves=4):
    """Generates MIDI numbers for a Dorian scale."""
    # Dorian scale intervals: Root, M2, m3, P4, P5, M6, m7
    # In semitones: 0, 2, 3, 5, 7, 9, 10
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10}
    scale_notes = set()
    root_note_val = start_midi % 12
    
    # Generate notes around the typical piano range
    for midi_val in range(21, 109):
        if (midi_val - root_note_val) % 12 in dorian_intervals:
            scale_notes.add(midi_val)
    return sorted(list(scale_notes))

def snap_to_scale(midi_seq, scale_notes):
    """Snaps each note in a sequence to the closest note in a given scale."""
    snapped_seq = []
    for note in midi_seq:
        # Find the closest note in the scale
        min_dist = float('inf')
        closest_note = -1
        
        # Find the two scale notes the current note is between
        for i in range(len(scale_notes) - 1):
            if scale_notes[i] <= note <= scale_notes[i+1]:
                dist1 = abs(note - scale_notes[i])
                dist2 = abs(note - scale_notes[i+1])
                # Tie-breaking rule: round up
                if dist1 < dist2:
                    closest_note = scale_notes[i]
                else: # Includes the tie case (dist1 == dist2)
                    closest_note = scale_notes[i+1]
                break
        # If the note was already in the scale, it will be found
        if note in scale_notes:
            closest_note = note

        snapped_seq.append(closest_note)
    return snapped_seq
    
# Initial sequence of notes
initial_notes = ['A4', 'C5', 'F4', 'E4', 'G4', 'C4', 'B4', 'D4']
midi_sequence = [note_to_midi(n) for n in initial_notes]

# --- 1. Transpose up Major Second (2 semitones) ---
op1_seq = [n + 2 for n in midi_sequence]

# --- 2. Invert Around E4 ---
axis_e4 = note_to_midi('E4')
op2_seq = [(2 * axis_e4) - n for n in op1_seq]

# --- 3. Retrograde the sequence ---
op3_seq = op2_seq[::-1]

# --- 4. Augment the intervals between consecutive notes by adding 3 semitones ---
op4_seq = [op3_seq[0]]
for i in range(1, len(op3_seq)):
    original_interval = op3_seq[i] - op3_seq[i-1]
    augmented_interval = original_interval + 3
    new_note = op4_seq[-1] + augmented_interval
    op4_seq.append(new_note)

# --- 5. Change to Dorian Mode Starting from D4 ---
d4_midi = note_to_midi('D4')
# First, transpose the sequence so its first note is D4
transpose_amount = d4_midi - op4_seq[0]
transposed_for_modal = [n + transpose_amount for n in op4_seq]
# Second, get the D Dorian scale notes
d_dorian_notes = get_dorian_scale_notes(d4_midi)
# Third, snap each note to the closest scale note
op5_seq = snap_to_scale(transposed_for_modal, d_dorian_notes)

# --- 6. Transpose down Minor Third (3 semitones) ---
op6_seq = [n - 3 for n in op5_seq]

# --- 7. Invert around F4 ---
axis_f4 = note_to_midi('F4')
op7_seq = [(2 * axis_f4) - n for n in op6_seq]

# --- 8. Transposed up one Octave (12 semitones) ---
final_midi_seq = [n + 12 for n in op7_seq]

# Convert final MIDI sequence back to note names
final_notes_seq = [midi_to_note(n) for n in final_midi_seq]

# Print the final result as a space-separated string
print(" ".join(final_notes_seq))
