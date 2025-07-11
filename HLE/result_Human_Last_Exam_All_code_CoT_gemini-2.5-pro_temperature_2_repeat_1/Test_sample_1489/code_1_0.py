import math

def solve_music_sequence():
    """
    Solves the musical transformation problem step-by-step.
    """

    # Helper functions to convert between note names and MIDI numbers.
    NOTE_NAMES_SHARP = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    NOTE_TO_VAL = {note: i for i, note in enumerate(NOTE_NAMES_SHARP)}

    def note_to_midi(note_str):
        """Converts a note string like 'A4' to a MIDI number."""
        note_str = note_str.upper()
        if len(note_str) > 2 and (note_str[1] == '#'):
            name = note_str[:2]
            octave = int(note_str[2:])
        else:
            name = note_str[0]
            octave = int(note_str[1:])
        
        val = NOTE_TO_VAL[name]
        return (octave + 1) * 12 + val

    def midi_to_note(midi_num):
        """Converts a MIDI number to a note string using sharps."""
        octave = midi_num // 12 - 1
        note_index = midi_num % 12
        note_name = NOTE_NAMES_SHARP[note_index]
        return f"{note_name}{octave}"

    # Step 0: Initial sequence
    initial_notes = "A4 C5 F4 E4 G4 C4 B4 D4".split()
    midi_sequence = [note_to_midi(n) for n in initial_notes]

    # Operation 1: Transpose up Major Second (+2 semitones)
    midi_sequence = [n + 2 for n in midi_sequence]

    # Operation 2: Invert Around E4 (MIDI 64)
    axis_e4 = note_to_midi('E4')
    midi_sequence = [(2 * axis_e4) - n for n in midi_sequence]

    # Operation 3: Retrograde
    midi_sequence.reverse()

    # Operation 4: Augment intervals by +3 semitones
    first_note = midi_sequence[0]
    intervals = [midi_sequence[i] - midi_sequence[i-1] for i in range(1, len(midi_sequence))]
    augmented_intervals = [i + 3 for i in intervals]
    new_sequence = [first_note]
    for interval in augmented_intervals:
        new_sequence.append(new_sequence[-1] + interval)
    midi_sequence = new_sequence

    # Operation 5: Change to Dorian Mode Starting from D4
    d_dorian_midi_mod_12 = {2, 4, 5, 7, 9, 11, 0} # D, E, F, G, A, B, C
    mapped_sequence = []
    for n in midi_sequence:
        mod_n = n % 12
        if mod_n in d_dorian_midi_mod_12:
            mapped_sequence.append(n)
        else:
            # Find closest scale tone, rounding up on tie
            dist_down, dist_up = 1, 1
            while (mod_n - dist_down) % 12 not in d_dorian_midi_mod_12: dist_down += 1
            while (mod_n + dist_up) % 12 not in d_dorian_midi_mod_12: dist_up += 1
            
            if dist_up < dist_down:
                mapped_sequence.append(n + dist_up)
            else: # dist_down < dist_up or it's a tie, we round up
                mapped_sequence.append(n - dist_down) if dist_down < dist_up else mapped_sequence.append(n + dist_up)

    midi_sequence = mapped_sequence

    # Operation 6: Transpose down Minor Third (-3 semitones)
    midi_sequence = [n - 3 for n in midi_sequence]

    # Operation 7: Invert around F4 (MIDI 65)
    axis_f4 = note_to_midi('F4')
    midi_sequence = [(2 * axis_f4) - n for n in midi_sequence]

    # Operation 8: Transposed up one Octave (+12 semitones)
    midi_sequence = [n + 12 for n in midi_sequence]

    # Final step: Convert MIDI back to note names and print
    final_note_names = [midi_to_note(n) for n in midi_sequence]
    print(" ".join(final_note_names))

solve_music_sequence()