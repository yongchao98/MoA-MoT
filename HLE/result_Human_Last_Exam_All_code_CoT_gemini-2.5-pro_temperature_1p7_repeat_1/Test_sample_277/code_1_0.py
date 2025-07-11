from fractions import Fraction

def solve_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    using note-to-note just intonation.
    """
    # 1. Define the melody and starting frequency
    # The melody represents the story of leaving ("G E E F D D C...")
    # and returning ("...G F E D C G C").
    # G4 is the starting G, G3 is the G an octave below.
    notes = [
        'G4', 'E4', 'E4', 'F4', 'D4', 'D4', 'C4', 'D4', 'E4', 'F4', 'G4', 'G4',
        'G4', 'F4', 'E4', 'D4', 'C4', 'G3', 'C4'
    ]
    start_freq = Fraction(392, 1)

    # 2. Define helper mappings for notes and intervals
    note_semitones = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}

    def get_pitch(note_str):
        """Converts a note string like 'G4' to a numerical semitone value."""
        note_name = note_str[0]
        octave = int(note_str[1])
        # We define C4 as pitch 24 for convenience
        return note_semitones[note_name] + 12 * (octave - 4) + 24

    # Ratios for just intervals, indexed by semitone difference
    # {semitone_diff: (Ratio, Interval Name)}
    interval_ratios = {
        0:   (Fraction(1, 1),   "Unison"),
        1:   (Fraction(16, 15), "Minor Second Up"),
        -1:  (Fraction(15, 16), "Minor Second Down"),
        2:   (Fraction(9, 8),   "Major Second Up"),
        -2:  (Fraction(8, 9),   "Major Second Down"),
        -3:  (Fraction(5, 6),   "Minor Third Down"), # M3 up is 5/4, down is 4/5, but G to E is m3 down
        5:   (Fraction(4, 3),   "Perfect Fourth Up"),
        -5:  (Fraction(3, 4),   "Perfect Fourth Down"),
    }

    # 3. Perform the step-by-step calculation
    current_freq = start_freq
    print(f"Starting note is {notes[0]} with frequency {current_freq.numerator}/{current_freq.denominator} Hz.\n")

    for i in range(1, len(notes)):
        prev_note_name = notes[i-1]
        current_note_name = notes[i]

        prev_pitch = get_pitch(prev_note_name)
        current_pitch = get_pitch(current_note_name)
        semitone_diff = current_pitch - prev_pitch

        ratio, interval_name = interval_ratios[semitone_diff]
        
        prev_freq = current_freq
        current_freq *= ratio

        print(f"Step {i}:")
        print(f"Note transition: {prev_note_name} -> {current_note_name} ({interval_name})")
        print(f"Frequency calculation: {prev_freq.numerator}/{prev_freq.denominator} Hz * {ratio.numerator}/{ratio.denominator} = {current_freq.numerator}/{current_freq.denominator} Hz")
        print("-" * 30)

    # 4. Final Result
    final_numerator = current_freq.numerator
    final_denominator = current_freq.denominator

    print(f"\nThe final frequency of the last note ({notes[-1]}) is {final_numerator}/{final_denominator} Hz.")
    print("\nThe result [a, b] is:")
    print([final_numerator, final_denominator])

solve_frequency()
<<<[62720, 243]>>>