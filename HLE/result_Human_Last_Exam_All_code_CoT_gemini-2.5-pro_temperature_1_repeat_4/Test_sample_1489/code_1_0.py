import math

def solve_music_sequence():
    """
    Solves the music theory problem by applying eight operations to a sequence of notes.
    """

    # Helper functions to convert between note names and MIDI numbers (C4=60)
    def note_name_to_midi(note_name):
        """Converts note name like A4 to standard MIDI number."""
        NOTE_MAP = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
        note_part = note_name[:-1]
        octave = int(note_name[-1])
        pitch = NOTE_MAP[note_part[0].upper()]
        if len(note_part) > 1:
            if note_part[1] == '#':
                pitch += 1
            elif note_part[1] == 'b':
                pitch -= 1
        return 12 * (octave + 1) + pitch

    def midi_to_note_name(midi_number):
        """Converts a MIDI number to a note name, using sharps for accidentals."""
        NOTE_NAMES_SHARP = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
        octave = (midi_number // 12) - 1
        note_index = midi_number % 12
        return f"{NOTE_NAMES_SHARP[note_index]}{octave}"

    # --- Define the 8 Operations ---

    # 1. Transpose
    def transpose(notes, semitones):
        return [n + semitones for n in notes]

    # 2. Invert
    def invert(notes, axis_midi):
        return [2 * axis_midi - n for n in notes]

    # 3. Retrograde
    def retrograde(notes):
        return notes[::-1]

    # 4. Augment Intervals
    def augment_intervals(notes, semitones_to_add):
        if not notes:
            return []
        new_notes = [notes[0]]
        for i in range(len(notes) - 1):
            original_interval = notes[i+1] - notes[i]
            augmented_interval = original_interval + semitones_to_add
            new_notes.append(new_notes[-1] + augmented_interval)
        return new_notes

    # 5. Change to Dorian Mode
    def snap_to_d_dorian(notes):
        # D Dorian scale notes are D, E, F, G, A, B, C
        # Pitch classes (C=0): 2, 4, 5, 7, 9, 11, 0
        d_dorian_pitch_classes = {0, 2, 4, 5, 7, 9, 11}
        snapped_notes = []
        for note in notes:
            offset = 0
            while True:
                upper_note = note + offset
                lower_note = note - offset
                is_upper_in_scale = (upper_note % 12) in d_dorian_pitch_classes
                is_lower_in_scale = (lower_note % 12) in d_dorian_pitch_classes

                # Tie-breaking rule: round up to the higher pitch
                if is_upper_in_scale and is_lower_in_scale and upper_note != lower_note:
                    snapped_notes.append(upper_note)
                    break
                elif is_upper_in_scale:
                    snapped_notes.append(upper_note)
                    break
                elif is_lower_in_scale:
                    snapped_notes.append(lower_note)
                    break
                offset += 1
        return snapped_notes

    # --- Main Execution ---

    # Initial sequence
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    current_midi = [note_name_to_midi(n) for n in initial_notes_str.split()]

    # Apply the 8 operations in order
    # 1. Transpose up Major Second (+2 semitones)
    current_midi = transpose(current_midi, 2)

    # 2. Invert Around E4 (MIDI 64)
    e4_midi = note_name_to_midi("E4")
    current_midi = invert(current_midi, e4_midi)

    # 3. Retrograde the sequence
    current_midi = retrograde(current_midi)

    # 4. Augment intervals by adding 3 semitones
    current_midi = augment_intervals(current_midi, 3)

    # 5. Change to Dorian Mode Starting from D4
    current_midi = snap_to_d_dorian(current_midi)

    # 6. Transpose down Minor Third (-3 semitones)
    current_midi = transpose(current_midi, -3)

    # 7. Invert around F4 (MIDI 65)
    f4_midi = note_name_to_midi("F4")
    current_midi = invert(current_midi, f4_midi)

    # 8. Transposed up one Octave (+12 semitones)
    current_midi = transpose(current_midi, 12)

    # Convert final MIDI numbers back to note names
    final_note_names = [midi_to_note_name(n) for n in current_midi]
    
    # Print the final result
    print(" ".join(final_note_names))

solve_music_sequence()