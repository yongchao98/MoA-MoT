from fractions import Fraction

def calculate_final_frequency():
    """
    This script calculates the frequency of the last note of "Hänschen klein".
    The calculation starts with a G at 392 Hz and applies just intonation ratios
    for each subsequent interval in the melody.
    """
    # 1. Define the melody and initial frequency.
    melody = ['G', 'E', 'E', 'F', 'D', 'D', 'C', 'D', 'E', 'F', 'G', 'G', 'G',
              'G', 'E', 'E', 'F', 'D', 'D', 'C', 'E', 'G', 'G', 'C', 'C', 'C']
    current_freq = Fraction(392, 1)

    # 2. Define the just intonation ratios for the intervals in the song.
    # The ratio is applied by multiplying it with the frequency of the preceding note.
    interval_ratios = {
        ('G', 'E'): Fraction(5, 6),    # Descending Minor Third
        ('E', 'E'): Fraction(1, 1),    # Unison
        ('E', 'F'): Fraction(16, 15),  # Ascending Minor Second
        ('F', 'D'): Fraction(5, 6),    # Descending Minor Third
        ('D', 'D'): Fraction(1, 1),    # Unison
        ('D', 'C'): Fraction(8, 9),    # Descending Major Second
        ('C', 'D'): Fraction(9, 8),    # Ascending Major Second
        ('D', 'E'): Fraction(9, 8),    # Ascending Major Second
        ('F', 'G'): Fraction(9, 8),    # Ascending Major Second
        ('G', 'G'): Fraction(1, 1),    # Unison
        ('C', 'E'): Fraction(5, 4),    # Ascending Major Third
        ('E', 'G'): Fraction(6, 5),    # Ascending Minor Third
        ('G', 'C'): Fraction(2, 3),    # Descending Perfect Fifth
        ('C', 'C'): Fraction(1, 1)     # Unison
    }

    # 3. Build the equation string while performing the calculation.
    # We start with the initial frequency.
    equation_str = str(current_freq.numerator)

    for i in range(len(melody) - 1):
        note1 = melody[i]
        note2 = melody[i+1]

        # Get the ratio for the interval.
        ratio = interval_ratios.get((note1, note2))

        # Update the frequency for the next note.
        current_freq *= ratio

        # Add the current multiplication to the equation string.
        # This fulfills the requirement to output each number in the final equation.
        equation_str += f" * ({ratio.numerator}/{ratio.denominator})"

    # 4. Get the final result as a coprime fraction a/b.
    final_a = current_freq.numerator
    final_b = current_freq.denominator

    # 5. Print the full equation and the final answer.
    print("The final frequency is the result of the following calculation:")
    print(f"{equation_str} = {final_a}/{final_b}")
    print("\nThe final answer as a list [a, b] is:")
    print([final_a, final_b])

calculate_final_frequency()