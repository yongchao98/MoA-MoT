from fractions import Fraction

def solve_song_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation.
    """
    # The starting note is G, tuned to 392 Hz.
    start_freq = Fraction(392, 1)

    # The melody for a verse of "Hänschen klein" is:
    # G E E | F D D | C D E F | G G G |
    # G E E | F D D | C E G G | C C C |
    #
    # The sequence of intervals and their corresponding just intonation ratios are:
    # G -> E: Down minor 3rd (5/6)
    # E -> E: Unison (1/1)
    # E -> F: Up minor 2nd (16/15)
    # F -> D: Down minor 3rd (5/6)
    # D -> D: Unison (1/1)
    # D -> C: Down major 2nd (8/9)
    # C -> D: Up major 2nd (9/8)
    # D -> E: Up major 2nd (9/8)
    # E -> F: Up minor 2nd (16/15)
    # F -> G: Up major 2nd (9/8)
    # G -> G: Unison (1/1) x3
    # G -> E: Down minor 3rd (5/6)
    # E -> E: Unison (1/1)
    # E -> F: Up minor 2nd (16/15)
    # F -> D: Down minor 3rd (5/6)
    # D -> D: Unison (1/1)
    # D -> C: Down major 2nd (8/9)
    # C -> E: Up major 3rd (5/4)
    # E -> G: Up minor 3rd (6/5)
    # G -> G: Unison (1/1)
    # G -> C: Down perfect 5th (2/3)
    # C -> C: Unison (1/1) x2

    ratios_str = [
        "5/6", "1/1", "16/15", "5/6", "1/1", "8/9", "9/8", "9/8", "16/15", "9/8",
        "1/1", "1/1", "1/1", "5/6", "1/1", "16/15", "5/6", "1/1", "8/9", "5/4",
        "6/5", "1/1", "2/3", "1/1", "1/1"
    ]

    ratios = [Fraction(s) for s in ratios_str]

    # Calculate the total ratio by multiplying all individual ratios
    total_ratio = Fraction(1, 1)
    for r in ratios:
        total_ratio *= r

    # The final frequency is the starting frequency multiplied by the total ratio.
    final_freq = start_freq * total_ratio

    a = final_freq.numerator
    b = final_freq.denominator

    # Output the components of the final equation
    print(f"Starting Frequency (G): {start_freq.numerator}/{start_freq.denominator} Hz")
    print(f"Total Frequency Ratio: {total_ratio.numerator}/{total_ratio.denominator}")
    print(f"Final Equation: ({start_freq.numerator}/{start_freq.denominator}) * ({total_ratio.numerator}/{total_ratio.denominator}) = {a}/{b}")
    
    # Print the final answer in the required list format
    print(f"The frequency of the last note is {a}/{b} Hertz.")
    print([a, b])

solve_song_frequency()