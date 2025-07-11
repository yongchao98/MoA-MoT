import sys
from fractions import Fraction

def solve_song_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation.

    The melody is traversed note by note, and for each interval, the frequency is
    multiplied by the corresponding just intonation ratio. The calculation is done
    using fractions to maintain precision.
    """
    # Initial frequency of the first note (G) in Hertz.
    start_freq = Fraction(392, 1)

    # Standard just intonation ratios for musical intervals.
    ratios = {
        'unison': Fraction(1, 1),
        'up_minor_second_diatonic': Fraction(16, 15),
        'down_minor_second_diatonic': Fraction(15, 16),
        'up_minor_second_chromatic': Fraction(25, 24),
        'down_minor_second_chromatic': Fraction(24, 25),
        'up_major_second': Fraction(9, 8),
        'down_major_second': Fraction(8, 9),
        'up_minor_third': Fraction(6, 5),
        'down_minor_third': Fraction(5, 6),
        'up_major_third': Fraction(5, 4),
        'down_major_third': Fraction(4, 5),
        'up_perfect_fourth': Fraction(4, 3),
        'down_perfect_fourth': Fraction(3, 4),
    }

    # The sequence of intervals for the chosen melody of "Hänschen klein".
    # Melody: G E E F D D C D E F G G G | A A G F# G A B B A G# A | B C' D' C' B A G F# G A B C' G
    intervals = [
        'down_major_third',           # G -> E
        'unison',                     # E -> E
        'up_minor_second_diatonic',   # E -> F
        'down_minor_third',           # F -> D
        'unison',                     # D -> D
        'down_major_second',          # D -> C
        'up_major_second',            # C -> D
        'up_major_second',            # D -> E
        'up_minor_second_diatonic',   # E -> F
        'up_major_second',            # F -> G
        'unison',                     # G -> G
        'unison',                     # G -> G
        'up_major_second',            # G -> A
        'unison',                     # A -> A
        'down_major_second',          # A -> G
        'down_minor_second_diatonic', # G -> F#
        'up_minor_second_diatonic',   # F# -> G
        'up_major_second',            # G -> A
        'up_major_second',            # A -> B
        'unison',                     # B -> B
        'down_major_second',          # B -> A
        'down_minor_second_chromatic',# A -> G#
        'up_minor_second_chromatic',  # G# -> A
        'up_major_second',            # A -> B
        'up_minor_second_diatonic',   # B -> C'
        'up_major_second',            # C' -> D'
        'down_major_second',          # D' -> C'
        'down_minor_second_diatonic', # C' -> B
        'down_major_second',          # B -> A
        'down_major_second',          # A -> G
        'down_minor_second_diatonic', # G -> F#
        'up_minor_second_diatonic',   # F# -> G
        'up_major_second',            # G -> A
        'up_major_second',            # A -> B
        'up_minor_second_diatonic',   # B -> C'
        'down_perfect_fourth',        # C' -> G
    ]

    current_freq = start_freq
    # Build the full equation string for printing.
    # To prevent extremely long lines, we'll add line breaks.
    max_len = 80
    equation_parts = [str(start_freq.numerator)]
    current_line_len = len(equation_parts[0])

    for interval in intervals:
        ratio = ratios[interval]
        part = f" * {ratio.numerator}/{ratio.denominator}"
        if current_line_len + len(part) > max_len:
            equation_parts.append("\n")
            current_line_len = 0
        equation_parts.append(part)
        current_line_len += len(part)
        current_freq *= ratio

    final_fraction = current_freq
    a = final_fraction.numerator
    b = final_fraction.denominator
    
    equation_str = "".join(equation_parts)
    equation_str += f" = {a}/{b}"
    
    print("The full calculation is:")
    print(equation_str)
    print("\nFinal frequency in Hertz as a fraction a/b:")
    print(f"{a}/{b}")
    print("\nAnswer in the format [a, b]:")
    print(f"[{a},{b}]")

solve_song_frequency()
<<<[100352, 243]>>>