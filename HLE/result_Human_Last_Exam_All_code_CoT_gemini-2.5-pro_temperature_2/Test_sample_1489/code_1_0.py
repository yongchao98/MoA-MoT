import re

def solve_music_sequence():
    """
    Solves the music theory problem by applying a series of transformations to a sequence of notes.
    """
    # --- Helper Functions ---

    NOTE_TO_MIDI_BASE = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    MIDI_TO_NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    def note_to_midi(note_str):
        """Converts a note string (e.g., 'A4', 'C5#') to a MIDI number."""
        match = re.match(r'([A-G])([#b]?)(\d+)', note_str)
        if not match:
            raise ValueError(f"Invalid note format: {note_str}")
        
        note_name, accidental, octave_str = match.groups()
        octave = int(octave_str)
        
        midi_val = NOTE_TO_MIDI_BASE[note_name]
        if accidental == '#':
            midi_val += 1
        elif accidental == 'b':
            midi_val -= 1
            
        return (octave + 1) * 12 + midi_val

    def midi_to_note(midi_num):
        """Converts a MIDI number to a note string (e.g., 'A4', 'C5#')."""
        if not isinstance(midi_num, int):
            raise TypeError("MIDI number must be an integer.")
            
        octave = (midi_num // 12) - 1
        note_name = MIDI_TO_NOTE_NAMES[midi_num % 12]
        
        if '#' in note_name:
            return f"{note_name[0]}{octave}{note_name[1]}"
        else:
            return f"{note_name}{octave}"

    def midi_sequence_to_notes(midi_seq):
        """Converts a sequence of MIDI numbers to a list of note strings."""
        return [midi_to_note(n) for n in midi_seq]

    # --- Initial Data ---
    initial_notes_str = "A4 C5 F4 E4 G4 C4 B4 D4"
    initial_notes_list = initial_notes_str.split()
    notes = [note_to_midi(n) for n in initial_notes_list]
    
    print(f"Initial Sequence: {initial_notes_list}")
    print(f"Initial MIDI: {notes}\n")

    # --- Operation 1: Transpose up Major Second ---
    notes = [n + 2 for n in notes]
    print("1. Transpose up Major Second (+2 semitones)")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")

    # --- Operation 2: Invert Around E4 ---
    pivot_e4 = note_to_midi('E4')
    notes = [pivot_e4 + (pivot_e4 - n) for n in notes]
    print(f"2. Invert Around E4 (MIDI {pivot_e4})")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")
    
    # --- Operation 3: Retrograde ---
    notes.reverse()
    print("3. Retrograde the sequence")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")
    
    # --- Operation 4: Augment Intervals by +3 semitones ---
    intervals = [notes[i] - notes[i-1] for i in range(1, len(notes))]
    augmented_intervals = [i + 3 for i in intervals]
    new_notes = [notes[0]]
    for interval in augmented_intervals:
        new_notes.append(new_notes[-1] + interval)
    notes = new_notes
    print("4. Augment intervals by adding 3 semitones")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")

    # --- Operation 5: Change to Dorian Mode Starting from D4 ---
    root_d4 = note_to_midi('D4')
    dorian_intervals = {0, 2, 3, 5, 7, 9, 10}
    dorian_pitch_classes = sorted([(root_d4 + i) % 12 for i in dorian_intervals])
    
    new_notes = []
    for note in notes:
        note_pc = note % 12
        if note_pc in dorian_pitch_classes:
            new_notes.append(note)
        else:
            # Find closest pitch class, rounding down on tie
            min_dist = float('inf')
            closest_pc = -1
            for scale_pc in dorian_pitch_classes:
                dist = abs(note_pc - scale_pc)
                dist = min(dist, 12 - dist)
                if dist < min_dist:
                    min_dist = dist
                    closest_pc = scale_pc
                elif dist == min_dist:
                    # Tie-break rule: round down, so prefer smaller pitch class
                    closest_pc = min(closest_pc, scale_pc)
            
            diff = closest_pc - note_pc
            if diff > 6: diff -= 12
            if diff < -6: diff += 12
            new_notes.append(note + diff)
    notes = new_notes
    print(f"5. Change to Dorian Mode Starting from D4 (Root MIDI {root_d4})")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")
    
    # --- Operation 6: Transpose down Minor Third ---
    notes = [n - 3 for n in notes]
    print("6. Transpose down Minor Third (-3 semitones)")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")
    
    # --- Operation 7: Invert around F4 ---
    pivot_f4 = note_to_midi('F4')
    notes = [pivot_f4 + (pivot_f4 - n) for n in notes]
    print(f"7. Invert around F4 (MIDI {pivot_f4})")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")
    
    # --- Operation 8: Transpose up one Octave ---
    notes = [n + 12 for n in notes]
    print("8. Transpose up one Octave (+12 semitones)")
    print(f"Resulting MIDI: {notes}")
    print(f"Resulting Notes: {' '.join(midi_sequence_to_notes(notes))}\n")
    
    # --- Final Answer ---
    final_note_sequence = ' '.join(midi_sequence_to_notes(notes))
    print("Final answer sequence:")
    print(final_note_sequence)
    print(f"<<<{final_note_sequence}>>>")

solve_music_sequence()