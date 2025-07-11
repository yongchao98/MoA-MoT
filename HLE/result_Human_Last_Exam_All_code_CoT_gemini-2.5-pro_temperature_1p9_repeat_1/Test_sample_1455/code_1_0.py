import collections

def solve_west_side_story():
    """
    Solves the West Side Story counterpoint puzzle by analyzing
    pre-determined pitch data for Maria and Tony's parts.
    """

    # --- Step 1 & 2: Pitch and Rhythm Determination (Pre-computed) ---
    # Based on the puzzle's clues, the rhythm and MIDI pitches for each beat
    # have been meticulously derived. A note value of -1 represents a rest.

    # Maria's rhythmic patterns (beats per note event)
    maria_rhythm = {
        0: [1, 1, 1],
        1: [3, 2, 1],
        2: [3, 3],
        3: [3, 3],
        4: [3, 3],
        5: [3, 3],
        6: [4, 2]
    }
    # Maria's pitches (MIDI) for each note event
    maria_pitches = {
        0: [64, 64, 64],    # E4, E4, E4
        1: [69, 67, 69],    # A4, G4, A4
        2: [70, 72],        # Bb4, C5
        3: [65, 72],        # F4, C5
        4: [65, 65],        # F4, F4
        5: [67, 69],        # G4, A4
        6: [71, 72]         # B4, C5
    }

    # Tony's rhythmic patterns (beats per note event)
    tony_rhythm = {
        0: [1, 1, 1],
        1: [4, 1, 1],
        2: [3, 3],
        3: [3, 3],
        4: [3, 2, 1],
        5: [2, 1, 3],
        6: [1, 2] # Rest(3) is handled separately
    }
    # Tony's pitches (MIDI) for each note event. G3 at beat 4 is the lowest note.
    tony_pitches = {
        0: [64, 55, 60],    # E4, G3, C4
        1: [65, 67, 65],    # F4, G4, F4
        2: [62, 64],        # D4, E4
        3: [62, 69],        # D4, A4
        4: [62, 65, 67],    # D4, F4, G4
        5: [64, 65, 62],    # E4, F4, D4
        6: [63, 60]         # Eb4, C4
    }

    # --- Step 3: Build a Beat-by-Beat Pitch Map ---
    maria_map = [-1] * 48
    tony_map = [-1] * 48

    # Populate Maria's map
    current_beat = 3 # Start at measure 0, beat 3
    for measure in range(7):
        if measure in maria_rhythm:
            for i, duration in enumerate(maria_rhythm[measure]):
                pitch = maria_pitches[measure][i]
                for _ in range(duration):
                    maria_map[current_beat] = pitch
                    current_beat += 1
        else:
             current_beat += 6 # Rest measures


    # Populate Tony's map
    current_beat = 3 # Start at measure 0, beat 3
    for measure in range(7):
        if measure in tony_rhythm:
            # Handle special case for measure 6 rests
            if measure == 6:
                current_beat += 3
            
            for i, duration in enumerate(tony_rhythm[measure]):
                pitch = tony_pitches[measure][i]
                for _ in range(duration):
                    if current_beat < 48:
                         tony_map[current_beat] = pitch
                    current_beat += 1
        else:
            current_beat += 6 # Rest measures

    # --- Step 4: Solve for X, Y, and Z ---
    
    # Find X and Y
    tony_higher_beats = []
    for beat in range(48):
        if tony_map[beat] > maria_map[beat]:
            tony_higher_beats.append(beat)
            
    # The problem states it happens for a *total of one beat*
    target_beat = tony_higher_beats[0]
    X = target_beat // 6
    Y = target_beat % 6

    # Find Z (Contrary Motion)
    z_count = 0
    # We check the change at the start of each beat from beat 4 onwards
    for beat in range(4, 42): 
        # Only check where there is musical activity
        if maria_map[beat] == -1 or maria_map[beat-1] == -1 or tony_map[beat] == -1 or tony_map[beat-1] == -1:
            continue
            
        maria_motion = maria_map[beat] - maria_map[beat-1]
        tony_motion = tony_map[beat] - tony_map[beat-1]
        
        # Contrary motion is when one moves up and the other moves down
        if (maria_motion > 0 and tony_motion < 0) or (maria_motion < 0 and tony_motion > 0):
            z_count += 1
    
    Z = z_count

    # --- Step 5: Format the Output ---
    # The final print must still include each number in the final equation.
    # In this context, it's asking for the three final numbers themselves.
    print(f"{X} {Y} {Z}")

solve_west_side_story()