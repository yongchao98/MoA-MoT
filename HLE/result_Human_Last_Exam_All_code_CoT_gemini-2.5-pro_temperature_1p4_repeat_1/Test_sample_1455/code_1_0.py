def solve_west_side_story():
    """
    Analyzes the counterpoint from West Side Story to find X, Y, and Z.
    """
    # Step 1 & 2: Decode notation and assign MIDI pitches.
    # MIDI Pitch Representation: C4=60, C#4=61, etc. None for rests.
    # Maria's part
    maria_pitches = [
        # m0: x x x E E E
        [None, None, None, 64, 64, 64],
        # m1: A---- G--
        [69, 69, 69, 69, 67, 67],
        # m2: B♭---- C--
        [70, 70, 70, 70, 72, 72],
        # m3: F---- C--
        [65, 65, 65, 65, 72, 72],
        # m4: F---- F--
        [65, 65, 65, 65, 65, 65],
        # m5: G---- A--
        [67, 67, 67, 67, 69, 69],
        # m6: B♮---- C-- (held)
        [71, 71, 71, 71, 72, 72],
        # m7: (held C)
        [72, 72, 72, 72, 72, 72],
    ]

    # Tony's part
    tony_pitches = [
        # m0: x x x E G C
        [None, None, None, 64, 55, 60], # G3 is lowest note, per clue
        # m1: F---- G F (aligned)
        [65, 65, 65, 65, 67, 65],
        # m2: D---- E--
        [62, 62, 62, 62, 64, 64],
        # m3: D---- A--
        [62, 62, 62, 62, 69, 69],
        # m4: D(4) F(1) G(1) interpretation for Tony > Maria on one beat
        [62, 62, 62, 62, 65, 67],
        # m5: E-- F D---
        [64, 64, 65, 62, 62, 62],
        # m6: --- E♭ C-- (held)
        [None, None, None, 63, 60, 60],
        # m7: (held C)
        [60, 60, 60, 60, 60, 60],
    ]

    # Flatten the lists for easier iteration
    m_flat = [beat for measure in maria_pitches for beat in measure]
    t_flat = [beat for measure in tony_pitches for beat in measure]

    # Step 3: Find X and Y
    # Find the single beat where Tony's voice is higher than Maria's.
    cross_over_beats = []
    for i in range(len(m_flat)):
        if m_flat[i] is not None and t_flat[i] is not None:
            if t_flat[i] > m_flat[i]:
                cross_over_beats.append(i)

    # Assuming there's only one such beat as per the prompt
    target_beat_index = cross_over_beats[0]
    X = target_beat_index // 6
    Y = target_beat_index % 6

    # Step 4: Find Z (Contrary Motion)
    # Contrary motion occurs when voices move in opposite directions.
    Z = 0
    for i in range(1, len(m_flat)):
        m_curr = m_flat[i]
        m_prev = m_flat[i - 1]
        t_curr = t_flat[i]
        t_prev = t_flat[i - 1]

        # Ensure we are comparing actual notes, not rests
        if all(p is not None for p in [m_curr, m_prev, t_curr, t_prev]):
            delta_m = m_curr - m_prev
            delta_t = t_curr - t_prev

            # Check for contrary motion (one goes up, other goes down)
            if (delta_m > 0 and delta_t < 0) or (delta_m < 0 and delta_t > 0):
                Z += 1

    # Step 5: Final Output
    print(f"{X} {Y} {Z}")

solve_west_side_story()
<<<4 5 4>>>