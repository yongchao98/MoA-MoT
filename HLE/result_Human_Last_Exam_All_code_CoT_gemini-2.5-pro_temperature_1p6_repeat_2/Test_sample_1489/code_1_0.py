import math

def solve_music_sequence():
    """
    Processes a sequence of musical notes through eight distinct operations.
    """

    # --- Helper Functions ---

    NOTE_PITCH_CLASSES = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    SHARP_NOTE_NAMES = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}

    def note_to_midi(note_str):
        """Converts a note name like 'A4' or 'C#5' to a MIDI number."""
        name = note_str[0].upper()
        octave_start_idx = 1
        accidental = 0
        if len(note_str) > 1:
            if note_str[1] == '#':
                accidental = 1
                octave_start_idx = 2
            elif note_str[1] == 'b':
                accidental = -1
                octave_start_idx = 2
        
        octave = int(note_str[octave_start_idx:])
        pitch_class = NOTE_PITCH_CLASSES[name] + accidental
        midi_val = 12 * (octave + 1) + pitch_class
        return midi_val

    def midi_to_note(midi_val):
        """Converts a MIDI number to a note name using sharps for accidentals."""
        if not isinstance(midi_val, int):
            raise TypeError("MIDI value must be an integer.")
        octave = (midi_val // 12) - 1
        pitch_class = midi_val % 12
        note_name = SHARP_NOTE_NAMES[pitch_class]
        return f"{note_name}{octave}"

    # --- Main Logic ---

    # Initial sequence
    initial_notes = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    current_midi = [note_to_midi(n) for n in initial_notes]

    # Operation 1: Transpose up Major Second (+2 semitones)
    current_midi = [n + 2 for n in current_midi]

    # Operation 2: Invert Around E4 (MIDI 64)
    pivot_e4 = note_to_midi('E4')
    current_midi = [(2 * pivot_e4) - n for n in current_midi]

    # Operation 3: Retrograde the sequence
    current_midi = current_midi[::-1]

    # Operation 4: Augment intervals by adding 3 semitones to each
    augmented_midi = []
    for i, note in enumerate(current_midi):
        # new_note_k = old_note_k + augmentation_amount * (k-1)
        augmented_midi.append(note + 3 * i)
    current_midi = augmented_midi

    # Operation 5: Change to Dorian Mode Starting from D4
    # D Dorian scale notes are D,E,F,G,A,B,C. Their pitch classes (mod 12) are:
    # C=0, D=2, E=4, F=5, G=7, A=9, B=11
    dorian_pcs = {0, 2, 4, 5, 7, 9, 11}
    
    def find_nearest_dorian_note(midi_note):
        pc = midi_note % 12
        if pc in dorian_pcs:
            return midi_note
        
        # Find distance to lower and upper scale notes
        dist_down = 1
        while (pc - dist_down) % 12 not in dorian_pcs:
            dist_down += 1
            
        dist_up = 1
        while (pc + dist_up) % 12 not in dorian_pcs:
            dist_up += 1

        # Return the closer note, rounding down in case of a tie
        if dist_down <= dist_up:
            return midi_note - dist_down
        else:
            return midi_note + dist_up

    current_midi = [find_nearest_dorian_note(n) for n in current_midi]
    
    # Operation 6: Transpose down Minor Third (-3 semitones)
    current_midi = [n - 3 for n in current_midi]

    # Operation 7: Invert around F4 (MIDI 65)
    pivot_f4 = note_to_midi('F4')
    current_midi = [(2 * pivot_f4) - n for n in current_midi]

    # Operation 8: Transpose up one Octave (+12 semitones)
    current_midi = [n + 12 for n in current_midi]

    # Convert final MIDI numbers back to note names and print
    final_notes = [midi_to_note(n) for n in current_midi]
    print(" ".join(final_notes))

if __name__ == "__main__":
    solve_music_sequence()