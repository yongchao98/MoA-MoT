from fractions import Fraction

def solve_song_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using
    sequential just intonation.
    """

    # 1. Deconstruct the Melody
    # The melody for one verse is broken into four standard phrases.
    # Note names with octave numbers ensure correct interval calculation.
    p1 = ['G4', 'E4', 'E4', 'F4', 'D4', 'D4', 'C4', 'D4', 'E4', 'F4', 'G4', 'G4', 'G4']
    p2 = ['G4', 'E4', 'E4', 'F4', 'D4', 'D4', 'C4', 'E4', 'G4', 'G4', 'C4']
    p3 = ['D4', 'D4', 'D4', 'D4', 'E4', 'F4', 'E4', 'E4', 'E4', 'E4', 'F4', 'G4']
    p4 = ['G4', 'E4', 'E4', 'F4', 'D4', 'D4', 'C4', 'E4', 'G4', 'G4', 'C4']
    
    # The full melody connects these phrases sequentially.
    # The interval between phrases is determined by the last note of the previous
    # phrase and the first note of the current one.
    melody_notes = p1 + p2 + p3 + p4

    # Helper function to convert note names to a numeric value (like MIDI)
    # to make interval calculation straightforward.
    note_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    def note_to_semitone(note_str):
        note_name = note_str[0]
        octave = int(note_str[1])
        # Using a standard reference (C0 = 12)
        return 12 * octave + note_map[note_name]

    melody_semitones = [note_to_semitone(n) for n in melody_notes]

    # 2. Define Just Intonation Ratios
    # Ratios for ascending intervals based on their size in semitones.
    just_ratios = {
        0: Fraction(1, 1),      # Unison
        1: Fraction(16, 15),    # Minor Second
        2: Fraction(9, 8),      # Major Second
        3: Fraction(6, 5),      # Minor Third
        4: Fraction(5, 4),      # Major Third
        5: Fraction(4, 3),      # Perfect Fourth
        7: Fraction(3, 2)       # Perfect Fifth
    }

    # 3. Calculate Cumulative Frequency Change
    # The starting note G is tuned to 392 Hz.
    initial_freq = Fraction(392, 1)
    current_freq = initial_freq
    
    print(f"Starting note: {melody_notes[0]}, Initial frequency: {initial_freq.numerator} Hz")
    
    # Calculate the product of all interval ratios.
    total_ratio = Fraction(1, 1)
    for i in range(1, len(melody_semitones)):
        prev_note_val = melody_semitones[i-1]
        curr_note_val = melody_semitones[i]
        
        interval = curr_note_val - prev_note_val
        
        ratio = Fraction(1, 1)
        if interval >= 0:  # Ascending interval or unison
            ratio = just_ratios[interval]
        else:  # Descending interval
            ratio = Fraction(1, just_ratios[abs(interval)])
        
        total_ratio *= ratio

    # 4. Final Calculation and Output
    final_freq = initial_freq * total_ratio

    print(f"\nInitial Frequency: {initial_freq.numerator}/{initial_freq.denominator}")
    print(f"Cumulative Ratio: {total_ratio.numerator}/{total_ratio.denominator}")
    print("\nFinal Frequency = Initial Frequency * Cumulative Ratio")
    print(f"Final Frequency = {initial_freq.numerator} * {total_ratio.numerator}/{total_ratio.denominator} = {final_freq.numerator}/{final_freq.denominator} Hz")
    
    final_answer = [final_freq.numerator, final_freq.denominator]
    print(f"\nThe final frequency of the last note ({melody_notes[-1]}) is {final_freq.numerator}/{final_freq.denominator} Hertz.")
    print(f"The result [a, b] is: {final_answer}")
    return final_answer

# Execute the function to get the answer.
result = solve_song_frequency()
print(f"\n<<<{result}>>>")