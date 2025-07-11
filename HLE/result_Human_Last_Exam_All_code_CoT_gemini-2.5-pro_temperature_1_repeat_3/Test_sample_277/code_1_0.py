from fractions import Fraction

def solve_song_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation
    for each interval relative to the preceding note.
    """

    # The note sequence for a full verse of "Hänschen klein" (Otto Frömmel version).
    # This sequence corresponds to the lyrics:
    # "Hänschen klein ging allein in die weite Welt hinein.
    # Stock und Hut steh'n ihm gut, ist gar wohlgemut.
    # Aber Mutter weinet sehr, hat ja nun kein Hänschen mehr!
    # Da besinnt sich das Kind, kehrt nach Haus geschwind."
    note_sequence = [
        'G4', 'E4', 'E4', 'F4', 'D4', 'D4', 'C4', 'D4', 'E4', 'F4', 'G4', 'G4', 'G4',
        'G4', 'E4', 'E4', 'F4', 'D4', 'D4', 'C4', 'E4', 'G4', 'G4', 'C4', 'C4', 'C4',
        'B3', 'C4', 'D4', 'E4', 'F4', 'F4', 'E4', 'D4', 'C4', 'B3', 'A3', 'G3', 'G3',
        'G3', 'B3', 'C4', 'D4', 'E4', 'F4', 'F4', 'E4', 'D4', 'G3', 'G3', 'C4', 'C4', 'C4'
    ]

    # Just intonation ratios for musical intervals.
    ratios = {
        "up m2": Fraction(16, 15), "down m2": Fraction(15, 16),
        "up M2": Fraction(9, 8),   "down M2": Fraction(8, 9),
        "up m3": Fraction(6, 5),   "down m3": Fraction(5, 6),
        "up M3": Fraction(5, 4),   "down M3": Fraction(4, 5),
        "up P4": Fraction(4, 3),   "down P4": Fraction(3, 4),
        "up P5": Fraction(3, 2),   "down P5": Fraction(2, 3),
        "unison": Fraction(1, 1)
    }

    # Helper maps for interval calculation.
    note_letters = ['C', 'D', 'E', 'F', 'G', 'A', 'B']
    letter_to_degree = {letter: i for i, letter in enumerate(note_letters)}
    letter_to_semitone = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}

    def get_interval_ratio(n1, n2):
        """
        Determines the interval between two notes and returns the just intonation ratio.
        """
        name1, oct1_str = n1[0], n1[1:]
        name2, oct2_str = n2[0], n2[1:]
        
        # Handle notes like 'G#4' if they were present
        oct1 = int(oct1_str.replace('#', '').replace('b', ''))
        oct2 = int(oct2_str.replace('#', '').replace('b', ''))

        semitone1 = 12 * oct1 + letter_to_semitone[name1]
        semitone2 = 12 * oct2 + letter_to_semitone[name2]
        s_diff = semitone2 - semitone1

        if s_diff == 0:
            return ratios["unison"]

        # Calculate degree difference from letter names
        deg1 = letter_to_degree[name1]
        deg2 = letter_to_degree[name2]
        
        # Determine interval number (e.g., 2nd, 3rd, 4th)
        if s_diff > 0: # Ascending interval
            interval_num = (deg2 - deg1 + 7) % 7 + 1
        else: # Descending interval
            interval_num = (deg1 - deg2 + 7) % 7 + 1
        
        direction = "up" if s_diff > 0 else "down"
        abs_s_diff = abs(s_diff)
        
        interval_name = ""
        if interval_num == 2: # Second
            interval_name = "M2" if abs_s_diff == 2 else "m2"
        elif interval_num == 3: # Third
            interval_name = "M3" if abs_s_diff == 4 else "m3"
        elif interval_num == 4: # Fourth
            # In just intonation, augmented 4th (tritone) can have other ratios, but doesn't occur here.
            interval_name = "P4"
        elif interval_num == 5: # Fifth
            interval_name = "P5"
        else: # Other intervals not present in this simple melody
            return Fraction(1, 1)

        key = f"{direction} {interval_name}"
        return ratios[key]

    # Starting frequency of G4 is 392 Hz
    start_freq = Fraction(392, 1)

    # Calculate total multiplier by chaining the ratios for each interval
    total_multiplier = Fraction(1, 1)
    for i in range(len(note_sequence) - 1):
        ratio = get_interval_ratio(note_sequence[i], note_sequence[i+1])
        total_multiplier *= ratio
    
    # Calculate final frequency using the start frequency and total multiplier
    final_freq = start_freq * total_multiplier

    # The final equation is: final_freq = start_freq * multiplier
    print(f"The equation for the final frequency is:")
    print(f"f_final = f_initial * M")
    print(f"where f_initial = {start_freq.numerator}/{start_freq.denominator} Hz")
    print(f"and the cumulative multiplier M = {total_multiplier.numerator}/{total_multiplier.denominator}")
    
    a = final_freq.numerator
    b = final_freq.denominator
    
    print(f"\nThe calculated frequency of the last note is:")
    print(f"{a} / {b} Hertz")
    
    print(f"\nThe final answer is the list [a,b]:")
    print(f"[{a},{b}]")
    print(f'<<<[{a},{b}]>>>')

solve_song_frequency()