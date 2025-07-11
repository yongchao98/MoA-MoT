import sys

def solve_counterpoint():
    """
    Analyzes the West Side Story counterpoint to find X, Y, and Z.
    """
    # Step 1: Establish Pitches (Note to MIDI mapping)
    # Based on the clues, the following octave assignments are determined.
    # R = Rest
    note_map = {
        'G3': 55, 'C4': 60, 'D4': 62, 'E♭4': 63, 'E4': 64, 'F4': 65,
        'G4': 67, 'A4': 69, 'B♭4': 70, 'B♮4': 71, 'C5': 72, 'R': None
    }

    # Step 2: Transcribe Melodies beat-by-beat
    # Each list represents the note sung on each of the 48 beats (8 measures * 6 beats)
    maria_notes = [
        'R', 'R', 'R', 'E4', 'E4', 'E4',  # Measure 0
        'A4', 'A4', 'A4', 'A4', 'G4', 'G4',  # Measure 1
        'B♭4', 'B♭4', 'B♭4', 'B♭4', 'C5', 'C5',  # Measure 2
        'F4', 'F4', 'F4', 'F4', 'C5', 'C5',  # Measure 3
        'F4', 'F4', 'F4', 'F4', 'F4', 'F4',  # Measure 4
        'G4', 'G4', 'G4', 'G4', 'A4', 'A4',  # Measure 5
        'B♮4', 'B♮4', 'B♮4', 'B♮4', 'C5', 'C5',  # Measure 6
        'R', 'R', 'R', 'R', 'R', 'R'   # Measure 7
    ]

    tony_notes = [
        'R', 'R', 'R', 'E4', 'G4', 'C4',  # Measure 0 (C is C4, not C5, to satisfy the "one beat higher" clue)
        'F4', 'F4', 'F4', 'F4', 'G4', 'F4',  # Measure 1
        'D4', 'D4', 'D4', 'D4', 'E4', 'E4',  # Measure 2
        'D4', 'D4', 'D4', 'D4', 'A4', 'A4',  # Measure 3
        'D4', 'D4', 'D4', 'D4', 'F4', 'G3',  # Measure 4 (G is G3 to satisfy the "lowest note" clue)
        'E4', 'F4', 'D4', 'D4', 'D4', 'D4',  # Measure 5
        'R', 'R', 'R', 'E♭4', 'C4', 'C4',  # Measure 6
        'R', 'R', 'R', 'R', 'R', 'R'   # Measure 7
    ]

    maria_midi = [note_map[n] for n in maria_notes]
    tony_midi = [note_map[n] for n in tony_notes]

    # Step 3: Find X and Y
    # Find the single beat where Tony's voice is higher than Maria's.
    x, y = -1, -1
    for i in range(len(maria_midi)):
        if maria_midi[i] is not None and tony_midi[i] is not None:
            if tony_midi[i] > maria_midi[i]:
                x = i // 6  # Measure number
                y = i % 6   # Beat number
                break

    # Step 4: Find Z
    # Count instances of contrary motion.
    z = 0
    # Iterate through note changes, not every beat.
    # A note change is defined as a point where a new note begins.
    note_change_indices = sorted(list(set(
        [6, 10, 12, 16, 18, 22, 24, 28, 30, 34, 36, 37, 38, 42, 45, 46]
    )))
    
    # We check the motion between points of melodic change.
    # A simpler way is to check every beat transition where both voices move.
    for i in range(1, len(maria_midi)):
        # Both voices must be singing at this beat and the previous one
        if maria_midi[i] is None or maria_midi[i-1] is None or \
           tony_midi[i] is None or tony_midi[i-1] is None:
            continue

        # Determine direction of movement for each voice
        maria_move = maria_midi[i] - maria_midi[i-1]
        tony_move = tony_midi[i] - tony_midi[i-1]

        # Both must move for it to be contrary motion
        if maria_move == 0 or tony_move == 0:
            continue

        # Check if they move in opposite directions
        if (maria_move > 0 and tony_move < 0) or \
           (maria_move < 0 and tony_move > 0):
            z += 1

    # Step 5: Output the result
    print(f"{x} {y} {z}")

solve_counterpoint()