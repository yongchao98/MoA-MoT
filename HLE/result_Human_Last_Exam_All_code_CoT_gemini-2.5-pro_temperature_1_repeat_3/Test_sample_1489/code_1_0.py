import re
import bisect

def note_to_midi(note_str):
    """Converts a note name (e.g., 'A4', 'C#5') to a MIDI number."""
    pitch_class_map = {'C': 0, 'C#': 1, 'D': 2, 'D#': 3, 'E': 4, 'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'A': 9, 'A#': 10, 'B': 11}
    match = re.match(r'([A-G]#?)([0-9])', note_str)
    name, octave = match.groups()
    
    # C4 is MIDI 60. Formula: pitch_class_val + (octave + 1) * 12
    return pitch_class_map[name] + (int(octave) + 1) * 12

def midi_to_note(midi_num):
    """Converts a MIDI number to a note name (e.g., 'A4', 'C5#')."""
    midi_pitch_map = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    pitch_class = midi_pitch_map[midi_num % 12]
    octave = midi_num // 12 - 1
    
    # Format as requested, e.g., C5# instead of C#5
    if '#' in pitch_class:
        return f"{pitch_class[0]}{octave}{pitch_class[1]}"
    else:
        return f"{pitch_class}{octave}"

def quantize_to_scale(notes, scale_tones):
    """Quantizes a list of notes to the nearest note in a given scale."""
    quantized_notes = []
    for note in notes:
        # Find insertion point to locate the two surrounding scale notes
        pos = bisect.bisect_left(scale_tones, note)
        
        # Handle edge cases
        if pos == 0:
            best_match = scale_tones[0]
        elif pos == len(scale_tones):
            best_match = scale_tones[-1]
        else:
            # Compare distances to the lower and upper scale notes
            low = scale_tones[pos - 1]
            high = scale_tones[pos]
            
            diff_low = abs(note - low)
            diff_high = abs(note - high)

            if diff_high < diff_low:
                best_match = high
            elif diff_low < diff_high:
                best_match = low
            else: # On tie, round up to the higher note
                best_match = high
        quantized_notes.append(best_match)
    return quantized_notes


def main():
    # Initial sequence of notes
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    notes = [note_to_midi(n) for n in initial_notes_str.split()]
    
    # 1. Transpose up Major Second (+2 semitones)
    notes = [n + 2 for n in notes]

    # 2. Invert Around E4 (MIDI 64)
    e4_midi = note_to_midi('E4')
    notes = [(2 * e4_midi) - n for n in notes]

    # 3. Retrograde the sequence
    notes = notes[::-1]
    
    # 4. Augment intervals by adding 3 semitones
    if len(notes) > 1:
        intervals = [notes[i+1] - notes[i] for i in range(len(notes) - 1)]
        augmented_intervals = [i + 3 for i in intervals]
        new_notes = [notes[0]]
        for interval in augmented_intervals:
            new_notes.append(new_notes[-1] + interval)
        notes = new_notes

    # 5. Change to Dorian Mode Starting from D4
    d4_midi = note_to_midi('D4')
    dorian_intervals = [0, 2, 3, 5, 7, 9, 10] # T-S-T-T-T-S-T in semitones from root
    # Generate D Dorian scale tones across a wide range
    d_dorian_scale_tones = []
    for octave in range(-1, 10):
        root = d4_midi + (octave - 4) * 12
        d_dorian_scale_tones.extend([root + i for i in dorian_intervals])
    notes = quantize_to_scale(notes, sorted(d_dorian_scale_tones))

    # 6. Transpose down Minor Third (-3 semitones)
    notes = [n - 3 for n in notes]

    # 7. Invert around F4 (MIDI 65)
    f4_midi = note_to_midi('F4')
    notes = [(2 * f4_midi) - n for n in notes]

    # 8. Transposed up one Octave (+12 semitones)
    notes = [n + 12 for n in notes]

    # Convert final MIDI numbers back to note names and print
    final_note_names = [midi_to_note(n) for n in notes]
    print(" ".join(final_note_names))

if __name__ == '__main__':
    main()
<<<A5 D6 C5# E5 B4 A4 C5# F4#>>>