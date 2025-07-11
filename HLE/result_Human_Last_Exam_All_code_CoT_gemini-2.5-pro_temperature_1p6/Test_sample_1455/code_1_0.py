def solve_west_side_story_counterpoint():
    """
    Analyzes the provided musical excerpt from West Side Story to find X, Y, and Z.

    - (X, Y): The measure and beat where Tony's voice is higher than Maria's.
    - Z: The total number of times the voices move in contrary motion.
    """

    # Step 1: Represent the Music as MIDI pitch lists.
    # None represents a rest. Each list represents 7 measures of 6 beats each.
    # This transcription is derived by interpreting the notation and applying all clues.
    # C4 = 60.
    
    # Maria: x x x E E E |A---- G-- |B♭---- C---- |F---- C---- |F---- F---- |G---- A---- |B♮---- C-----|
    # Rhythmic interpretation (beats per note):
    # M0: rest(3),1,1,1 | M1: 4,2 | M2: 3,3 | M3: 3,3 | M4: 6 | M5: 4,2 | M6: 4,2
    maria_pitches = [
        # M0         # M1              # M2              # M3
        None,None,None,64,64,64, 69,69,69,69,67,67, 70,70,70,72,72,72, 65,65,65,72,72,72,
        # M4              # M5              # M6
        65,65,65,65,65,65, 67,67,67,67,69,69, 71,71,71,71,72,72
    ]

    # Tony: x x x E G C |F------ G F |D---- E---- |D---- A---- |D---- F-- G |E-- F D-----|--- E♭ C-----|
    # Rhythmic interpretation (beats per note):
    # M0: rest(3),1,1,1 | M1: 4,1,1 | M2: 3,3 | M3: 3,3 | M4: 3,2,1 | M5: 2,1,3 | M6: rest(3),1,2
    tony_pitches = [
        # M0         # M1              # M2              # M3
        None,None,None,64,67,60, 65,65,65,65,67,65, 62,62,62,64,64,64, 62,62,62,69,69,69,
        # M4              # M5              # M6
        62,62,62,65,65,55, 64,64,65,62,62,62, None,None,None,63,60,60
    ]

    # --- Find X and Y ---
    X, Y = -1, -1
    for i in range(len(maria_pitches)):
        if tony_pitches[i] is not None and maria_pitches[i] is not None:
            if tony_pitches[i] > maria_pitches[i]:
                X = i // 6
                Y = i % 6
                # As per the clue, this happens only once.
                break

    # --- Find Z ---
    Z = 0
    # We start checking from the second beat of actual music
    # (index 4, which is beat 4 of measure 0)
    for i in range(4, len(maria_pitches)):
        # Get pitches at the current and previous beat
        m_curr, m_prev = maria_pitches[i], maria_pitches[i-1]
        t_curr, t_prev = tony_pitches[i], tony_pitches[i-1]
        
        # Skip if there's a rest or the start of music
        if None in [m_curr, m_prev, t_curr, t_prev]:
            continue

        # Calculate the change in pitch (motion) for each voice
        maria_motion = m_curr - m_prev
        tony_motion = t_curr - t_prev

        # Contrary motion is when one motion is positive and the other is negative
        if maria_motion * tony_motion < 0:
            Z += 1
            
    # Print the final result in the format "X Y Z"
    print(f"{X} {Y} {Z}")

solve_west_side_story_counterpoint()
<<<0 4 4>>>