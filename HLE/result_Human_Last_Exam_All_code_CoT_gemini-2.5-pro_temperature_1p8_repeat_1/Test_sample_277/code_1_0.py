import math
from fractions import Fraction

def solve_hanschen_klein_frequency():
    """
    Calculates the final frequency of the song "Hänschen klein"
    using just intonation intervals applied sequentially for each note in the melody.
    """
    
    # The melody of "Hänschen klein" (Otto Frömmel version).
    # 'H' represents the high C (c'') to distinguish it from the standard C (c').
    melody_str = "G E E F D D C D E F G G G G E E F D D C E G G F D C D E F G A G F E D C D E F G F E D G A B H B A G F E D C"
    melody = melody_str.split()

    # The starting frequency of G is 392 Hz. We use fractions for precision.
    start_freq = Fraction(392, 1)
    
    print("This script calculates the frequency of the final note of 'Hänschen klein'.")
    print(f"The song starts on note {melody[0]} with a frequency of {start_freq.numerator}/{start_freq.denominator} Hz.")
    print("The frequency of each subsequent note is calculated from the previous one using just intonation ratios.")
    print("-" * 30)

    # Map note names to semitone values relative to C=0.
    note_map_semitones = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11, 'H': 12}
    
    current_freq = start_freq
    print(f"Start: Note {melody[0]}, Frequency = {start_freq.numerator}/{start_freq.denominator} Hz")
    
    # Iterate through the melody to calculate the frequency of each note
    for i in range(len(melody) - 1):
        n1_str = melody[i]
        n2_str = melody[i+1]
        
        n1_val = note_map_semitones[n1_str]
        n2_val = note_map_semitones[n2_str]
        
        semitone_diff = n2_val - n1_val

        # Map semitone difference to just intonation frequency ratios.
        ratio = Fraction(1, 1)
        interval_name = "Unknown"
        if semitone_diff == 0: 
            ratio, interval_name = Fraction(1, 1), "Unison"
        elif semitone_diff == 1: 
            ratio, interval_name = Fraction(16, 15), "Asc. Minor 2nd"
        elif semitone_diff == -1: 
            ratio, interval_name = Fraction(15, 16), "Desc. Minor 2nd"
        elif semitone_diff == 2: 
            ratio, interval_name = Fraction(9, 8), "Asc. Major 2nd"
        elif semitone_diff == -2: 
            ratio, interval_name = Fraction(8, 9), "Desc. Major 2nd"
        elif semitone_diff == 3: 
            ratio, interval_name = Fraction(6, 5), "Asc. Minor 3rd"
        elif semitone_diff == -3: 
            ratio, interval_name = Fraction(5, 6), "Desc. Minor 3rd"
        elif semitone_diff == 4: 
            ratio, interval_name = Fraction(5, 4), "Asc. Major 3rd"
        elif semitone_diff == -4: 
            ratio, interval_name = Fraction(4, 5), "Desc. Major 3rd"
        elif semitone_diff == 5: 
            ratio, interval_name = Fraction(4, 3), "Asc. Perfect 4th"
        elif semitone_diff == -5: 
            ratio, interval_name = Fraction(3, 4), "Desc. Perfect 4th"

        prev_freq_num, prev_freq_den = current_freq.numerator, current_freq.denominator
        current_freq *= ratio
        
        print(
            f"Step {i+1:2d}: {n1_str} -> {n2_str} ({interval_name}). Ratio = {ratio.numerator}/{ratio.denominator}. "
            f"Equation: {prev_freq_num}/{prev_freq_den} * {ratio.numerator}/{ratio.denominator} = {current_freq.numerator}/{current_freq.denominator} Hz"
        )
        
    print("-" * 30)
    final_freq = current_freq
    print(f"The frequency of the last note ({melody[-1]}) is {final_freq.numerator}/{final_freq.denominator} Hertz.")
    
    # The Fraction class automatically keeps the fraction in simplest form (coprime).
    final_answer = [final_freq.numerator, final_freq.denominator]
    print(f"\nThe final answer in the format [a, b] is: {final_answer}")
    
    # Returning the final answer as requested for extraction.
    return final_answer

# Execute the function to get the answer
final_result = solve_hanschen_klein_frequency()
print(f'<<<{final_result}>>>')
