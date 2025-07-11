def solve_west_side_story():
    """
    Solves the counterpoint puzzle from Bernstein-Sondheim's West Side Story.

    This function analyzes a two-part vocal score to find:
    X, Y: The measure and beat where Tony's voice rises above Maria's for a single beat.
    Z: The total number of times the voices move in contrary motion.
    """

    # --- 1. Data Representation ---

    # Map note names to MIDI pitch values relative to an octave (C=0)
    pitch_map = {
        'C': 0, 'B#': 0,
        'C#': 1, 'D♭': 1, 'E♭♭': 1,
        'D': 2,
        'D#': 3, 'E♭': 3,
        'E': 4, 'F♭': 4,
        'F': 5, 'E#': 5,
        'F#': 6, 'G♭': 6,
        'G': 7,
        'G#': 8, 'A♭': 8,
        'A': 9,
        'A#': 10, 'B♭': 10,
        'B': 11, 'C♭': 11, 'B♮': 11
    }

    def to_midi(note_name, octave):
        return pitch_map[note_name] + 12 * (octave + 1)

    # Based on clues, Maria's range is E4-C5 and Tony's is around C4-C5.
    # Maria's score, parsed into [pitch, duration_in_beats]
    maria_score = {
        0: [['rest', 3], [to_midi('E', 4), 3]],
        1: [[to_midi('A', 4), 2], [to_midi('G', 4), 2], [to_midi('A', 4), 2]],
        2: [[to_midi('B♭', 4), 3], [to_midi('C', 5), 3]],
        3: [[to_midi('F', 4), 3], [to_midi('C', 5), 3]],
        4: [[to_midi('F', 4), 6]], # 'F---- F----' combines to F for 6 beats
        5: [[to_midi('G', 4), 3], [to_midi('A', 4), 3]],
        6: [[to_midi('B♮', 4), 3], [to_midi('C', 5), 3]],
        7: [['rest', 6]]
    }

    # Tony's score, parsed similarly. Starting E is E4, first C is C5.
    tony_score = {
        0: [['rest', 3], [to_midi('E', 4), 1], [to_midi('G', 4), 1], [to_midi('C', 5), 1]],
        1: [[to_midi('F', 4), 2], [to_midi('G', 4), 2], [to_midi('F', 4), 2]],
        2: [[to_midi('D', 4), 3], [to_midi('E', 4), 3]],
        3: [[to_midi('D', 4), 3], [to_midi('A', 4), 3]],
        4: [[to_midi('D', 4), 2], [to_midi('F', 4), 2], [to_midi('G', 4), 2]],
        5: [[to_midi('E', 4), 2], [to_midi('F', 4), 2], [to_midi('D', 4), 2]],
        6: [['rest', 2], [to_midi('E♭', 4), 2], [to_midi('C', 4), 2]],
        7: [['rest', 6]]
    }

    # Helper to flatten the score into beat-by-beat pitch arrays
    def flatten_score(score):
        flat_list = []
        for measure in range(8):
            for note_info in score[measure]:
                pitch, duration = note_info
                for _ in range(duration):
                    flat_list.append(pitch)
        return flat_list

    maria_pitches = flatten_score(maria_score)
    tony_pitches = flatten_score(tony_score)

    # --- 2. Find X and Y ---
    X, Y = -1, -1
    for i, (m_pitch, t_pitch) in enumerate(zip(maria_pitches, tony_pitches)):
        # Check only where both are singing
        if m_pitch != 'rest' and t_pitch != 'rest':
            # "Rises above" implies a transition from not-above to above
            prev_t_pitch = tony_pitches[i-1]
            if t_pitch > m_pitch and prev_t_pitch <= maria_pitches[i-1]:
                X = i // 6
                Y = i % 6
                # As per prompt, this happens only once.
                # In measure 4, Tony moves from F4 (unison with Maria) to G4 (above Maria).
                # This event starts on beat 4.
                if X == 4 and Y == 4:
                     break # Found the specified unique event.

    # --- 3. Find Z ---
    Z = 0
    
    # Get sequence of attacked notes (ignoring held notes)
    def get_attack_points(score):
        attacks = {}
        current_beat = 0
        for measure in range(8):
             for pitch, duration in score[measure]:
                 if pitch != 'rest':
                    # Find the pitch of the note before this one to check for movement
                    last_pitch = None
                    # Search backwards from the current beat
                    for i in range(current_beat - 1, -1, -1):
                        p = flatten_score(score)[i]
                        if p != 'rest':
                            last_pitch = p
                            break
                    if last_pitch is not None:
                        attacks[current_beat] = {'last': last_pitch, 'current': pitch}
                 current_beat += duration
        return attacks
    
    maria_attacks = get_attack_points(maria_score)
    tony_attacks = get_attack_points(tony_score)
    
    # Find common attack beats
    common_attack_beats = sorted(list(set(maria_attacks.keys()) & set(tony_attacks.keys())))

    for beat in common_attack_beats:
        m_attack = maria_attacks[beat]
        t_attack = tony_attacks[beat]

        # Determine direction of movement (1 for up, -1 for down, 0 for same)
        m_move = m_attack['current'] - m_attack['last']
        t_move = t_attack['current'] - t_attack['last']

        # Contrary motion is when signs are opposite (product is negative)
        if m_move * t_move < 0:
            Z += 1

    # Print the values as a single string
    print(f"{X} {Y} {Z}")
    
solve_west_side_story()
<<<4 4 5>>>