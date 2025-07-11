from fractions import Fraction

def solve_song_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    using just intonation from the preceding note.
    """
    # 1. Starting frequency and note sequence
    start_freq = 392  # G4 is 392 Hz
    # Melody based on the "geschwind" version
    # G E E | F D D | C D E F | G G G ||
    # G E E | F D D | C E G G ||
    notes = [
        'G', 'E', 'E', 'F', 'D', 'D', 'C', 'D', 'E', 'F', 'G', 'G', 'G',
        'G', 'E', 'E', 'F', 'D', 'D', 'C', 'E', 'G', 'G'
    ]

    # 2. Define just intonation ratios for the intervals present in the song
    ratios = {
        # Unison
        ('G', 'G'): Fraction(1, 1),
        ('E', 'E'): Fraction(1, 1),
        ('D', 'D'): Fraction(1, 1),
        ('C', 'C'): Fraction(1, 1),
        # Seconds
        ('E', 'F'): Fraction(16, 15), # Ascending Minor Second
        ('C', 'D'): Fraction(9, 8),  # Ascending Major Second
        ('D', 'E'): Fraction(9, 8),  # Ascending Major Second
        ('F', 'G'): Fraction(9, 8),  # Ascending Major Second
        ('D', 'C'): Fraction(8, 9),  # Descending Major Second
        # Thirds
        ('E', 'G'): Fraction(6, 5),  # Ascending Minor Third
        ('C', 'E'): Fraction(5, 4),  # Ascending Major Third
        ('G', 'E'): Fraction(4, 5),  # Descending Major Third
        ('F', 'D'): Fraction(4, 5),  # Descending Major Third
    }

    # 3. Calculate the frequency by chaining the ratios
    current_freq = Fraction(start_freq)
    
    # Build a string to display the calculation
    calculation_steps = [str(start_freq)]
    
    # Iterate through the note transitions
    for i in range(len(notes) - 1):
        note1 = notes[i]
        note2 = notes[i+1]
        
        transition = (note1, note2)
        if transition not in ratios:
            print(f"Error: Ratio for transition {transition} not found.")
            return

        ratio = ratios[transition]
        current_freq *= ratio
        
        # Add the operation to the calculation string
        if ratio.denominator == 1:
            calculation_steps.append(f"* {ratio.numerator}")
        else:
            calculation_steps.append(f"* ({ratio.numerator}/{ratio.denominator})")

    # 4. Print the final results
    print("The final frequency is calculated as follows:")
    final_equation = ' '.join(calculation_steps)
    
    # Let's print each number/fraction in the final equation as requested
    print("Final Equation:", end=" ")
    for step in calculation_steps:
      print(step, end=" ")
    print() # Newline

    print(f"\nFinal Result: {current_freq.numerator} / {current_freq.denominator} Hertz")
    
    # The final answer in the specified list format
    final_list = [current_freq.numerator, current_freq.denominator]
    print(final_list)


solve_song_frequency()

# The problem is solved, so here is the final answer in the required format.
# a = 2**19 * 7**2 = 25690112
# b = 5**7 = 78125
final_answer = [25690112, 78125]
print(f"<<<{final_answer}>>>")