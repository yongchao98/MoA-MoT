import math

def main():
    """
    This script performs a series of eight musical operations on a sequence of notes.
    """

    # --- Helper functions for note conversion ---
    NOTE_NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]

    def note_to_midi(note_str):
        """Converts a note string like 'A4' to a MIDI number."""
        note_str = note_str.strip().upper()
        name = note_str[:-1]
        octave = int(note_str[-1])
        
        pitch_class = NOTE_NAMES.index(name)
        return pitch_class + (octave + 1) * 12

    def midi_to_note(midi_num):
        """Converts a MIDI number to a note string like 'A4'."""
        octave = (midi_num // 12) - 1
        pitch_class = midi_num % 12
        note_name = NOTE_NAMES[pitch_class]
        return f"{note_name}{octave}"

    # --- Initial State ---
    # Initial sequence of musical notes: A4 C5 F4 E4 G4 C4 B4 D4
    # Corresponding MIDI values: 69 72 65 64 67 60 71 62
    current_sequence = [69, 72, 65, 64, 67, 60, 71, 62]

    # --- Apply the eight operations sequentially ---

    # 1. Transpose up Major Second (+2 semitones)
    current_sequence = [note + 2 for note in current_sequence]
    
    # 2. Invert Around E4 (MIDI 64)
    axis_e4 = 64
    current_sequence = [(2 * axis_e4) - note for note in current_sequence]

    # 3. Retrograde the sequence
    current_sequence.reverse()

    # 4. Augment the intervals by adding 3 semitones to each interval.
    if len(current_sequence) > 1:
        augmented_sequence = [current_sequence[0]]
        for i in range(len(current_sequence) - 1):
            original_interval = current_sequence[i+1] - current_sequence[i]
            new_interval = original_interval + 3
            next_note = augmented_sequence[-1] + new_interval
            augmented_sequence.append(next_note)
        current_sequence = augmented_sequence

    # 5. Change to Dorian Mode Starting from D4 (Root MIDI 62)
    # D Dorian scale pitch classes (relative to C=0): {D, E, F, G, A, B, C} -> {2, 4, 5, 7, 9, 11, 0}
    dorian_pitch_classes = {0, 2, 4, 5, 7, 9, 11}
    quantized_sequence = []
    for note in current_sequence:
        if note % 12 in dorian_pitch_classes:
            quantized_sequence.append(note)
        else:
            # Find the nearest scale note, tie-breaking by choosing the lower note
            offset = 1
            while True:
                lower_candidate = note - offset
                if lower_candidate % 12 in dorian_pitch_classes:
                    quantized_sequence.append(lower_candidate)
                    break
                upper_candidate = note + offset
                if upper_candidate % 12 in dorian_pitch_classes:
                    quantized_sequence.append(upper_candidate)
                    break
                offset += 1
    current_sequence = quantized_sequence

    # 6. Transpose down Minor Third (-3 semitones)
    current_sequence = [note - 3 for note in current_sequence]

    # 7. Invert around F4 (MIDI 65)
    axis_f4 = 65
    current_sequence = [(2 * axis_f4) - note for note in current_sequence]

    # 8. Transposed up one Octave (+12 semitones)
    current_sequence = [note + 12 for note in current_sequence]

    # --- Final Output ---
    # Convert final MIDI sequence back to note names for the final answer.
    final_note_names = [midi_to_note(n) for n in current_sequence]
    
    # Print the final sequence of eight notes as a space-separated string.
    print(" ".join(final_note_names))

if __name__ == "__main__":
    main()