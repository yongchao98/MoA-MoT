from fractions import Fraction

def calculate_final_frequency():
    """
    Calculates the frequency of the last note of a verse from "Hänschen klein"
    using just intonation based on the preceding note.
    """

    # The starting note is G, with a frequency of 392 Hz.
    # We represent this as a fraction for precise calculations.
    initial_freq = Fraction(392, 1)

    # The melody for the verse ending "geschwind" is G-E-E-F-D-D-C-E-G-G-C.
    # We list the ratios for the intervals between consecutive notes.
    # G -> E: Descending Major Third (4/5)
    # E -> E: Unison (1/1)
    # E -> F: Ascending Minor Second (16/15)
    # F -> D: Descending Major Third (4/5)
    # D -> D: Unison (1/1)
    # D -> C: Descending Major Second (8/9)
    # C -> E: Ascending Major Third (5/4)
    # E -> G: Ascending Minor Third (6/5)
    # G -> G: Unison (1/1)
    # G -> C: Descending Perfect Fifth (2/3)
    
    ratios = [
        Fraction(4, 5),    # G -> E
        Fraction(1, 1),    # E -> E
        Fraction(16, 15),  # E -> F
        Fraction(4, 5),    # F -> D
        Fraction(1, 1),    # D -> D
        Fraction(8, 9),    # D -> C
        Fraction(5, 4),    # C -> E
        Fraction(6, 5),    # E -> G
        Fraction(1, 1),    # G -> G
        Fraction(2, 3)     # G -> C (the last note)
    ]

    # Start with the frequency of the first note of the verse (G).
    current_freq = initial_freq

    # Multiply the frequency by each ratio in sequence.
    for ratio in ratios:
        current_freq *= ratio
    
    # The final frequency as a fraction of coprime integers a/b.
    a = current_freq.numerator
    b = current_freq.denominator

    # As per the instruction "output each number in the final equation".
    print(f"The final frequency is the result of the equation: 392 * (product of all ratios)")
    print(f"Final frequency = {a}/{b} Hz")
    
    # Print the answer in the specified list format [a, b].
    print([a, b])

calculate_final_frequency()