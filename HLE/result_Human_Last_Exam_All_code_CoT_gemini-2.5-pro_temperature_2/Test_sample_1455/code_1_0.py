def solve_counterpoint():
    """
    Solves the West Side Story counterpoint problem.
    """
    # Step 1 & 2: Parse score and assign pitches based on musical deduction.
    # The musical notation has been interpreted and pitches assigned according to the clues.
    # N (None) represents a rest. Pitches are standard MIDI note numbers.
    # Maria's highest note C5=72. Tony's lowest G3=55. Start in unison. End an octave apart.
    
    # Maria's pitches beat-by-beat for measures 0-7
    m_pitches = [
        # M0: x x x E E E
        None, None, None, 64, 64, 64,
        # M1: A---- G-- A
        69, 69, 69, 67, 67, 69,
        # M2: B♭---- C----
        70, 70, 70, 72, 72, 72,
        # M3: F---- C----
        65, 65, 65, 72, 72, 72,
        # M4: F---- F----
        65, 65, 65, 65, 65, 65,
        # M5: G---- A----
        67, 67, 67, 69, 69, 69,
        # M6: B♮---- C-----
        71, 71, 71, 72, 72, 72,
        # M7: ------------
        None, None, None, None, None, None
    ]

    # Tony's pitches beat-by-beat for measures 0-7
    t_pitches = [
        # M0: x x x E G C
        None, None, None, 64, 55, 60,
        # M1: F------ G F  (Interpreted as F(4 beats), G(1), F(1))
        65, 65, 65, 65, 67, 65,
        # M2: D---- E----
        62, 62, 62, 64, 64, 64,
        # M3: D---- A----
        62, 62, 62, 69, 69, 69,
        # M4: D---- F-- G
        62, 62, 62, 65, 65, 67,
        # M5: E-- F D-----
        64, 64, 65, 62, 62, 62,
        # M6: --- E♭ C-----
        None, None, 63, 60, 60, 60,
        # M7: ------------
        None, None, None, None, None, None
    ]

    # Step 3: Find X and Y
    # Find the single beat where Tony's voice is higher than Maria's.
    X, Y = -1, -1
    for i in range(len(m_pitches)):
        if m_pitches[i] is not None and t_pitches[i] is not None:
            if t_pitches[i] > m_pitches[i]:
                X = i // 6  # measure number
                Y = i % 6   # beat number
                break

    # Step 4: Find Z
    # Count occurrences of contrary motion.
    # Contrary motion happens when one voice moves up and the other moves down.
    Z = 0
    # We track the last heard pitch to calculate melodic movement.
    last_m_pitch = None
    last_t_pitch = None
    
    # We look for points where both melodies change pitch simultaneously.
    for i in range(len(m_pitches)):
        current_m_pitch = m_pitches[i]
        current_t_pitch = t_pitches[i]

        if current_m_pitch is not None and current_t_pitch is not None:
            # Check for motion only if there was a previous note to move from.
            if last_m_pitch is not None and last_t_pitch is not None:
                m_moved = current_m_pitch != last_m_pitch
                t_moved = current_t_pitch != last_t_pitch

                if m_moved and t_moved:
                    m_direction = current_m_pitch - last_m_pitch
                    t_direction = current_t_pitch - last_t_pitch
                    
                    # Contrary motion is when their products are negative.
                    if m_direction * t_direction < 0:
                        Z += 1

            # Update last known pitches
            last_m_pitch = current_m_pitch
            last_t_pitch = current_t_pitch

    # Step 5: Format and print the final answer
    print(f"{X} {Y} {Z}")

solve_counterpoint()