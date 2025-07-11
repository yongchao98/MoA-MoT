import math

def solve_west_side_story():
    """
    Solves the West Side Story counterpoint puzzle by parsing musical notation,
    assigning MIDI pitches, and analyzing the melodic relationship.
    """
    # Note to MIDI base value (octave added later)
    note_map = {
        'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11,
        'B♭': 10, 'E♭': 3, 'B♮': 11
    }

    # Assign octave (C4 = Middle C = 60)
    # Pitch assignment is based on clues: start in unison, end octave apart,
    # stay within an octave, high/low note limits, and minimal leaps.
    # A None value represents a rest.
    maria_midi = [
        None, None, None, 64, 64, 64,  # M0: x x x E E E (E4)
        69, 69, 69, 69, 67, 67,        # M1: A---- G-- (A4, G4)
        70, 70, 70, 72, 72, 72,        # M2: B♭-- C-- (B-flat4, C5)
        65, 65, 65, 72, 72, 72,        # M3: F-- C-- (F4, C5)
        65, 65, 65, 65, 65, 65,        # M4: F------ (F4)
        67, 67, 67, 69, 69, 69,        # M5: G-- A-- (G4, A4)
        71, 71, 71, 72, 72, 72,        # M6: B♮-- C--- (B4, C5)
        72, 72, 72, 72, 72, 72         # M7: (held C5)
    ]

    tony_midi = [
        None, None, None, 64, 67, 60,  # M0: x x x E G C (E4, G4, C4)
        65, 65, 65, 65, 65, 65,        # M1: F------ (F4)
        62, 62, 62, 64, 64, 64,        # M2: D-- E-- (D4, E4)
        62, 62, 62, 69, 69, 69,        # M3: D-- A-- (D4, A4)
        62, 62, 62, 62, 65, 55,        # M4: D---- F- G- (D4, F4, G3)
        64, 64, 65, 62, 62, 62,        # M5: E-- F- D--- (E4, F4, D4)
        None, None, None, 63, 60, 60,  # M6: --- E♭ C-- (E-flat4, C4)
        60, 60, 60, 60, 60, 60         # M7: (held C4)
    ]

    # --- Find X and Y (Tony's voice higher) ---
    x = -1
    y = -1
    higher_count = 0
    for i in range(len(maria_midi)):
        if maria_midi[i] is not None and tony_midi[i] is not None:
            if tony_midi[i] > maria_midi[i]:
                higher_count += 1
                x = i // 6
                y = i % 6

    # The problem states it happens for a total of one beat. Our analysis confirms this.

    # --- Find Z (Contrary Motion) ---
    z = 0
    for i in range(1, len(maria_midi)):
        m_curr, m_prev = maria_midi[i], maria_midi[i-1]
        t_curr, t_prev = tony_midi[i], tony_midi[i-1]

        # Skip if either voice is resting at the start or end of the interval
        if m_curr is None or m_prev is None or t_curr is None or t_prev is None:
            continue

        # Check if both voices move
        m_moved = m_curr != m_prev
        t_moved = t_curr != t_prev

        if m_moved and t_moved:
            # Determine direction of movement (positive for up, negative for down)
            m_dir = m_curr - m_prev
            t_dir = t_curr - t_prev

            # Check for contrary motion (opposite signs)
            if (m_dir > 0 and t_dir < 0) or (m_dir < 0 and t_dir > 0):
                z += 1
    
    print(f"{x} {y} {z}")

solve_west_side_story()
<<<0 4 2>>>