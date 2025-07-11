import math

def note_to_midi(note_str):
    """Converts a note name string (e.g., 'A4', 'C#5') to a MIDI number."""
    NOTES = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    
    # Find the split point between note name and octave number
    name = note_str.rstrip('0123456789')
    octave = int(note_str[len(name):])
    
    pitch_class = NOTES[name]
    return (octave + 1) * 12 + pitch_class

def midi_to_note(midi):
    """Converts a MIDI number to a note name string."""
    NOTES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    octave = midi // 12 - 1
    note_index = midi % 12
    return f"{NOTES[note_index]}{octave}"

def map_to_scale(notes, root, scale_intervals):
    """Maps each note to the nearest note in a given scale, choosing lower in case of a tie."""
    
    def is_in_scale(midi_note):
        return (midi_note - root) % 12 in scale_intervals

    new_notes = []
    for note in notes:
        if is_in_scale(note):
            new_notes.append(note)
            continue
        
        # Search outwards for the closest scale note
        d = 1
        while True:
            # Check lower note first to handle ties by choosing the lower pitch
            lower_note = note - d
            if is_in_scale(lower_note):
                new_notes.append(lower_note)
                break
            
            upper_note = note + d
            if is_in_scale(upper_note):
                new_notes.append(upper_note)
                break
            d += 1
    return new_notes

def main():
    # Initial sequence
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    midi_sequence = [note_to_midi(n) for n in initial_notes_str.split()]

    # 1. Transpose up Major Second (+2 semitones)
    seq1 = [n + 2 for n in midi_sequence]

    # 2. Invert Around E4 (MIDI 64)
    axis2 = note_to_midi('E4')
    seq2 = [axis2 + (axis2 - n) for n in seq1]

    # 3. Retrograde
    seq3 = list(reversed(seq2))

    # 4. Augment intervals (+3 semitones)
    seq4 = [seq3[0]]
    for i in range(1, len(seq3)):
        interval = seq3[i] - seq3[i-1]
        new_note = seq4[-1] + interval + 3
        seq4.append(new_note)

    # 5. Change to Dorian Mode Starting from D4
    root5 = note_to_midi('D4')
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10}  # Semitone steps from root: W-H-W-W-W-H-W
    seq5 = map_to_scale(seq4, root5, dorian_intervals)

    # 6. Transpose down Minor Third (-3 semitones)
    seq6 = [n - 3 for n in seq5]

    # 7. Invert around F4 (MIDI 65)
    axis7 = note_to_midi('F4')
    seq7 = [axis7 + (axis7 - n) for n in seq6]

    # 8. Transposed up one Octave (+12 semitones)
    final_midi_sequence = [n + 12 for n in seq7]

    # Convert final MIDI list to note names and print
    final_note_names = [midi_to_note(n) for n in final_midi_sequence]
    print(" ".join(final_note_names))

if __name__ == "__main__":
    main()