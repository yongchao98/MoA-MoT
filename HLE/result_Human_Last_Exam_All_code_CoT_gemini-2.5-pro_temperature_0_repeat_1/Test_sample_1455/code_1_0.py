def solve_counterpoint():
    """
    Analyzes the West Side Story counterpoint to find X, Y, and Z.
    X, Y: The measure and beat where Tony's voice is higher than Maria's.
    Z: The number of times contrary motion occurs.
    """

    # Step 1 & 2: The provided musical lines have been analyzed and converted to
    # MIDI pitch numbers based on the clues. A rest is represented by None.
    # C4 (Middle C) = 60.

    # Maria's pitches per beat (Measures 0-7)
    maria_pitches = [
        None, None, None, 64, 64, 64,  # M0: E4
        69, 69, 69, 69, 67, 67,  # M1: A4, G4
        70, 70, 70, 70, 72, 72,  # M2: Bb4, C5
        65, 65, 65, 65, 72, 72,  # M3: F4, C5
        65, 65, 65, 65, 65, 65,  # M4: F4
        67, 67, 67, 67, 69, 69,  # M5: G4, A4
        71, 71, 71, 71, 72, 72,  # M6: B4, C5
        None, None, None, None, None, None  # M7
    ]

    # Tony's pitches per beat (Measures 0-7), deduced from clues and notation
    tony_pitches = [
        None, None, None, 64, 55, 60,  # M0: E4, G3, C4
        65, 65, 67, 67, 65, 65,  # M1: F4, G4, F4
        62, 62, 62, 62, 64, 64,  # M2: D4, E4
        62, 62, 62, 62, 69, 69,  # M3: D4, A4
        62, 62, 62, 62, 65, 67,  # M4: D4, F4, G4
        64, 64, 65, 62, 62, 62,  # M5: E4, F4, D4
        None, None, None, 63, 60, 60,  # M6: Eb4, C4
        None, None, None, None, None, None  # M7
    ]

    # Step 3: Find X and Y
    # Iterate through each beat to find where Tony's pitch is higher than Maria's.
    X, Y = -1, -1
    for i in range(len(maria_pitches)):
        # Ensure both are singing (not a rest)
        if tony_pitches[i] is not None and maria_pitches[i] is not None:
            if tony_pitches[i] > maria_pitches[i]:
                X = i // 6  # Measure number
                Y = i % 6   # Beat number
                # The problem states this happens only once, so we can stop.
                break

    # Step 4: Find Z
    # Count instances of contrary motion.
    Z = 0
    for i in range(1, len(maria_pitches)):
        # Get current and previous pitches for both parts
        m_curr = maria_pitches[i]
        m_prev = maria_pitches[i-1]
        t_curr = tony_pitches[i]
        t_prev = tony_pitches[i-1]

        # Contrary motion requires both voices to change pitch from the previous beat.
        # Skip if there's a rest or if one of the voices holds a note.
        if m_curr is None or m_prev is None or t_curr is None or t_prev is None:
            continue
        if m_curr == m_prev or t_curr == t_prev:
            continue

        # Determine direction of movement for each voice
        maria_movement = m_curr - m_prev
        tony_movement = t_curr - t_prev

        # Check if they moved in opposite directions (one positive, one negative)
        if (maria_movement > 0 and tony_movement < 0) or \
           (maria_movement < 0 and tony_movement > 0):
            Z += 1

    # Step 5: Format and print the final answer
    print(f"{X} {Y} {Z}")

solve_counterpoint()
<<<4 5 3>>>