import fractions

def solve_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    using just intonation relative to the preceding note.
    """
    
    # 1. Set the initial frequency for G4.
    # We use the fractions.Fraction type to maintain precision.
    start_freq = fractions.Fraction(392)

    # 2. Define the melody for the first verse.
    melody = "G E E F D D C D E F G G G".split()

    # 3. Define the just intonation ratios for the intervals found in the melody.
    # The song is in C major, which has two sizes of major second:
    # Major Tone (9/8) and Minor Tone (10/9).
    ratios = {
        ('G', 'E'): fractions.Fraction(5, 6),    # Descending minor third
        ('E', 'E'): fractions.Fraction(1, 1),    # Unison
        ('E', 'F'): fractions.Fraction(16, 15),   # Ascending minor second (semitone)
        ('F', 'D'): fractions.Fraction(5, 6),    # Descending minor third
        ('D', 'D'): fractions.Fraction(1, 1),    # Unison
        ('D', 'C'): fractions.Fraction(8, 9),    # Descending major second (Major Tone)
        ('C', 'D'): fractions.Fraction(9, 8),    # Ascending major second (Major Tone)
        ('D', 'E'): fractions.Fraction(10, 9),   # Ascending major second (Minor Tone)
        # ('E', 'F') is already defined above
        ('F', 'G'): fractions.Fraction(9, 8),    # Ascending major second (Major Tone)
        ('G', 'G'): fractions.Fraction(1, 1)     # Unison
    }

    # 4. Sequentially calculate the frequency through the melody.
    current_freq = start_freq
    equation_parts = [str(start_freq.numerator)]

    for i in range(len(melody) - 1):
        note1 = melody[i]
        note2 = melody[i+1]
        
        # Look up the ratio for the current interval.
        try:
            ratio = ratios[(note1, note2)]
        except KeyError:
            print(f"Error: Ratio for interval ({note1}, {note2}) not defined.")
            return

        # Apply the ratio to the current frequency.
        current_freq *= ratio
        
        # Build the equation string for explanation.
        equation_parts.append(f"* ({ratio.numerator}/{ratio.denominator})")

    # 5. The final frequency is the result. The fraction is already simplified.
    final_freq = current_freq
    a = final_freq.numerator
    b = final_freq.denominator
    
    # 6. Print the results as requested.
    print("The final frequency is a/b Hertz, where a and b are coprime.")
    print("The full calculation is:")
    print(f"{' '.join(equation_parts)} = {a}/{b}")
    print("-" * 20)
    print("The numbers in the final equation are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print("-" * 20)
    print("The final answer in list format [a, b] is:")
    print([a, b])
    
    # Final answer in the specified format
    print(f"\n<<<[{a},{b}]>>>")

solve_frequency()