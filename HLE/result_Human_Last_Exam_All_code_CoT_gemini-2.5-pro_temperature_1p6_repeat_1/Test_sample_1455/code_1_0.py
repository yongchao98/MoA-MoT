import collections

def solve_counterpoint():
    """
    Analyzes a counterpoint from West Side Story to find specific musical events.
    """
    # Step 1 & 2: Define pitches and rhythms based on analysis
    # MIDI Pitch mapping. C4 is middle C.
    pitch_map = {
        'C3': 48, 'D3': 50, 'E3': 52, 'F3': 53, 'G3': 55, 'A3': 57, 'B3': 59,
        'C4': 60, 'D♭4': 61, 'D4': 62, 'E♭4': 63, 'E4': 64, 'F4': 65, 'G4': 67, 'A4': 69, 'B♭4': 70, 'B4': 71,
        'C5': 72,
    }
    
    # Rhythms and notes are interpreted from the puzzle's notation and clues.
    # Each inner list is a measure of 6 beats. A value of 0 represents a rest.
    maria_notes = [
        [0, 0, 0, pitch_map['E4'], pitch_map['E4'], pitch_map['E4']],
        [pitch_map['A4'], pitch_map['A4'], pitch_map['A4'], pitch_map['A4'], pitch_map['G4'], pitch_map['G4']],
        [pitch_map['B♭4'], pitch_map['B♭4'], pitch_map['B♭4'], pitch_map['C5'], pitch_map['C5'], pitch_map['C5']],
        [pitch_map['F4'], pitch_map['F4'], pitch_map['F4'], pitch_map['C5'], pitch_map['C5'], pitch_map['C5']], # C5 to avoid voice crossing
        [pitch_map['F4'], pitch_map['F4'], pitch_map['F4'], pitch_map['F4'], pitch_map['F4'], pitch_map['F4']],
        [pitch_map['G4'], pitch_map['G4'], pitch_map['G4'], pitch_map['A4'], pitch_map['A4'], pitch_map['A4']],
        [pitch_map['B4'], pitch_map['B4'], pitch_map['B4'], pitch_map['C5'], pitch_map['C5'], pitch_map['C5']],
        [pitch_map['C5'], pitch_map['C5'], pitch_map['C5'], pitch_map['C5'], pitch_map['C5'], pitch_map['C5']],
    ]

    tony_notes = [
        [0, 0, 0, pitch_map['E4'], pitch_map['G3'], pitch_map['C4']], # G3 satisfies "lowest note G" clue
        [pitch_map['F4'], pitch_map['F4'], pitch_map['F4'], pitch_map['F4'], pitch_map['G4'], pitch_map['F4']],
        [pitch_map['D4'], pitch_map['D4'], pitch_map['D4'], pitch_map['E4'], pitch_map['E4'], pitch_map['E4']],
        [pitch_map['D4'], pitch_map['D4'], pitch_map['D4'], pitch_map['A4'], pitch_map['A4'], pitch_map['A4']],
        [pitch_map['D4'], pitch_map['D4'], pitch_map['D4'], pitch_map['F4'], pitch_map['F4'], pitch_map['G4']],
        [pitch_map['E4'], pitch_map['E4'], pitch_map['F4'], pitch_map['D4'], pitch_map['D4'], pitch_map['D4']],
        [pitch_map['D4'], pitch_map['D4'], pitch_map['D4'], pitch_map['E♭4'], pitch_map['E♭4'], pitch_map['C4']],
        [pitch_map['C4'], pitch_map['C4'], pitch_map['C4'], pitch_map['C4'], pitch_map['C4'], pitch_map['C4']],
    ]
    
    # Flatten the lists of measures into a single sequence of beats
    maria_pitches = [pitch for measure in maria_notes for pitch in measure]
    tony_pitches = [pitch for measure in tony_notes for pitch in measure]
    
    # Step 3: Find X and Y
    # Find the single beat where Tony's voice is higher than Maria's.
    X, Y = -1, -1
    crossings = []
    for i in range(len(tony_pitches)):
        # Ignore rests
        if tony_pitches[i] > 0 and maria_pitches[i] > 0:
            if tony_pitches[i] > maria_pitches[i]:
                measure_num = i // 6
                beat_num = i % 6
                crossings.append((measure_num, beat_num))

    # The problem states there is exactly one beat of crossing.
    if len(crossings) == 1:
        X, Y = crossings[0]
        
    # Step 4: Find Z (Contrary Motion)
    # Count instances where one voice moves up and the other moves down.
    Z = 0
    for i in range(1, len(tony_pitches)):
        # Ensure we are not looking at a rest from the previous beat
        if tony_pitches[i-1] == 0 or maria_pitches[i-1] == 0:
            continue
            
        # Check if both voices change pitch
        is_tony_move = tony_pitches[i] != tony_pitches[i-1]
        is_maria_move = maria_pitches[i] != maria_pitches[i-1]

        if is_tony_move and is_maria_move:
            delta_tony = tony_pitches[i] - tony_pitches[i-1]
            delta_maria = maria_pitches[i] - maria_pitches[i-1]
            
            # Contrary motion if one delta is positive and the other is negative
            if delta_tony * delta_maria < 0:
                Z += 1

    # Step 5: Format and print the final answer
    print(f"{X} {Y} {Z}")

solve_counterpoint()
<<<4 5 4>>>