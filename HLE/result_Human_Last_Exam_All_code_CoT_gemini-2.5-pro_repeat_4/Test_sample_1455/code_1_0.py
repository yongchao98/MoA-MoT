import math

def solve_music_puzzle():
    """
    Solves the West Side Story counterpoint puzzle.
    """

    # 1. Define raw musical data
    maria_measures_str = [
        "x x x E E E",
        "A---- G-- A",
        "B♭---- C----",
        "F---- C----",
        "F---- F----",
        "G---- A----",
        "B♮---- C-----",
        "------------"
    ]
    tony_measures_str = [
        "x x x E G C",
        "F------ G F",
        "D---- E----",
        "D---- A----",
        "D---- F-- G",
        "E-- F D-----",
        "--- E♭ C-----",
        "------------"
    ]

    # 2. Define pitch mappings based on the clues
    # Clues analysis:
    # - Unison on E -> Maria's E and Tony's E are the same pitch. Let's start with E4 (MIDI 64).
    # - Maria's highest is C -> This must be C5 (MIDI 72).
    # - Tony's lowest is G -> This must be G3 (MIDI 55).
    # - End an octave apart -> Maria's final C is C5(72), Tony's is C4(60). Difference is 12 (one octave).
    # - Always within an octave -> This set of pitches holds this true.
    maria_pitch_map = {'E': 64, 'F': 65, 'G': 67, 'A': 69, 'B♭': 70, 'B♮': 71, 'C': 72}
    tony_pitch_map = {'C': 60, 'D': 62, 'E♭': 63, 'E': 64, 'F': 65, 'G': 55, 'A': 69}

    def parse_measure(measure_str, pitch_map):
        """
        Parses a measure string into a list of 6 pitches.
        """
        notes = measure_str.split()
        num_notes = len(notes)
        beats_per_note = 6 // num_notes
        
        pitch_list = []
        for note in notes:
            # Handle rests
            if note == 'x' or note.startswith('-'):
                pitch = None
            else:
                pitch = pitch_map.get(note)

            pitch_list.extend([pitch] * beats_per_note)
        return pitch_list

    # 3. Generate full pitch sequences
    maria_pitches = []
    for m_str in maria_measures_str:
        maria_pitches.extend(parse_measure(m_str, maria_pitch_map))

    tony_pitches = []
    for t_str in tony_measures_str:
        tony_pitches.extend(parse_measure(t_str, tony_pitch_map))

    # 4. Find X and Y (Tony's voice is higher)
    x = -1
    y = -1
    for i in range(len(maria_pitches)):
        if tony_pitches[i] is not None and maria_pitches[i] is not None:
            if tony_pitches[i] > maria_pitches[i]:
                x = i // 6
                y = i % 6
                break # Problem states this happens for only one beat

    # 5. Find Z (Contrary motion)
    z = 0
    for i in range(1, len(maria_pitches)):
        # Get pitches for the current and previous beat
        m_curr, m_prev = maria_pitches[i], maria_pitches[i-1]
        t_curr, t_prev = tony_pitches[i], tony_pitches[i-1]

        # A melodic motion requires two consecutive sounded notes
        if m_curr is None or m_prev is None or t_curr is None or t_prev is None:
            continue

        # Calculate direction of movement
        m_dir = m_curr - m_prev
        t_dir = t_curr - t_prev

        # Check for contrary motion (one moves up, the other down)
        if (m_dir > 0 and t_dir < 0) or (m_dir < 0 and t_dir > 0):
            z += 1
            
    # Print the final result
    print(f"{x} {y} {z}")

solve_music_puzzle()