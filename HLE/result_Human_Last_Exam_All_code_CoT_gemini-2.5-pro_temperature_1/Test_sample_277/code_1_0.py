from fractions import Fraction

def solve_haenschen_klein_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation
    for each interval relative to the preceding note.
    """
    # 1. Model the melody for one verse of the song.
    melody_notes = ['G', 'E', 'E', 'F', 'D', 'D', 'C', 'D', 'E', 'F', 'G', 'G', 'G',
                    'G', 'E', 'E', 'F', 'D', 'D', 'C', 'E', 'G', 'G', 'C', 'C', 'C']

    # 2. Define just intonation frequency ratios for the required intervals.
    # These ratios are derived from the C-major just scale (C=1, D=9/8, E=5/4, F=4/3, G=3/2).
    # The ratio for an interval from Note1 to Note2 is Freq(Note2)/Freq(Note1).
    interval_ratios = {
        ('G', 'E'): Fraction(5, 4) / Fraction(3, 2),   # Downward Minor 3rd -> 5/6
        ('E', 'F'): Fraction(4, 3) / Fraction(5, 4),   # Upward Minor 2nd -> 16/15
        ('F', 'D'): Fraction(9, 8) / Fraction(4, 3),   # Downward Minor 3rd -> 27/32
        ('D', 'C'): Fraction(1, 1) / Fraction(9, 8),   # Downward Major 2nd -> 8/9
        ('C', 'D'): Fraction(9, 8) / Fraction(1, 1),   # Upward Major 2nd -> 9/8
        ('D', 'E'): Fraction(5, 4) / Fraction(9, 8),   # Upward Major 2nd (minor tone) -> 10/9
        ('F', 'G'): Fraction(3, 2) / Fraction(4, 3),   # Upward Major 2nd -> 9/8
        ('C', 'E'): Fraction(5, 4) / Fraction(1, 1),   # Upward Major 3rd -> 5/4
        ('E', 'G'): Fraction(3, 2) / Fraction(5, 4),   # Upward Minor 3rd -> 6/5
        ('G', 'C'): Fraction(1, 1) / Fraction(3, 2),   # Downward Perfect 5th -> 2/3
    }

    # 3. Set the initial frequency of G as a Fraction.
    initial_freq = Fraction(392)
    current_freq = initial_freq
    
    # Store the ratios used for building the equation string
    used_ratios = []

    # 4. Iterate through the melody, calculating the frequency of each note.
    for i in range(len(melody_notes) - 1):
        start_note = melody_notes[i]
        end_note = melody_notes[i+1]
        
        # Unison interval has a ratio of 1/1
        if start_note == end_note:
            ratio = Fraction(1, 1)
        else:
            interval = (start_note, end_note)
            ratio = interval_ratios[interval]
        
        used_ratios.append(ratio)
        current_freq *= ratio

    final_freq = current_freq
    
    # 5. Construct and print the equation.
    equation_str = f"{initial_freq.numerator}"
    for r in used_ratios:
        # We only show non-unison intervals to keep the equation readable
        if r != 1:
            equation_str += f" * ({r.numerator}/{r.denominator})"

    # 6. Get the final coprime integers a and b.
    a = final_freq.numerator
    b = final_freq.denominator
    
    equation_str += f" = {a}/{b}"
    
    print("Calculation of the final frequency:")
    print(equation_str)
    
    # 7. Print the final answer in the required list format.
    print("\nThe final frequency is a/b Hertz, where [a, b] is:")
    print([a, b])

solve_haenschen_klein_frequency()