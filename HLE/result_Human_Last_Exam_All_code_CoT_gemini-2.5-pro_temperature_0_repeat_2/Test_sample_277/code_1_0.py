from fractions import Fraction

def calculate_final_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    based on a sequence of just intervals.
    """
    # The starting frequency of the first note (G) is 392 Hz.
    start_freq = Fraction(392, 1)

    # The sequence of notes for the verse is:
    # G, E, E, F, D, D, C, D, E, F, G, G, G,
    # G, E, E, F, D, D, C, E, G, G, C, C, C
    
    # The corresponding sequence of just interval ratios is determined.
    # For example, G down to E is a minor third (ratio 5/6),
    # E up to F is a minor second (ratio 16/15), etc.
    ratios = [
        Fraction(5, 6),    # G -> E: minor 3rd down
        Fraction(1, 1),    # E -> E: unison
        Fraction(16, 15),  # E -> F: minor 2nd up
        Fraction(5, 6),    # F -> D: minor 3rd down
        Fraction(1, 1),    # D -> D: unison
        Fraction(8, 9),    # D -> C: major 2nd down
        Fraction(9, 8),    # C -> D: major 2nd up
        Fraction(9, 8),    # D -> E: major 2nd up
        Fraction(16, 15),  # E -> F: minor 2nd up
        Fraction(9, 8),    # F -> G: major 2nd up
        Fraction(1, 1),    # G -> G: unison
        Fraction(1, 1),    # G -> G: unison
        Fraction(1, 1),    # G -> G: unison
        Fraction(5, 6),    # G -> E: minor 3rd down
        Fraction(1, 1),    # E -> E: unison
        Fraction(16, 15),  # E -> F: minor 2nd up
        Fraction(5, 6),    # F -> D: minor 3rd down
        Fraction(1, 1),    # D -> D: unison
        Fraction(8, 9),    # D -> C: major 2nd down
        Fraction(5, 4),    # C -> E: major 3rd up
        Fraction(6, 5),    # E -> G: minor 3rd up
        Fraction(1, 1),    # G -> G: unison
        Fraction(2, 3),    # G -> C: perfect 5th down
        Fraction(1, 1),    # C -> C: unison
        Fraction(1, 1)     # C -> C: unison
    ]

    # Calculate the total cumulative ratio by multiplying all individual ratios.
    total_ratio = Fraction(1, 1)
    for r in ratios:
        total_ratio *= r

    # Calculate the final frequency by applying the total ratio to the start frequency.
    final_freq = start_freq * total_ratio

    # Extract the numerator (a) and denominator (b) for the final answer.
    a = final_freq.numerator
    b = final_freq.denominator

    # Print the steps of the calculation and the final result.
    print("Calculation of the final note's frequency:")
    print(f"Initial frequency (G): {start_freq.numerator}/{start_freq.denominator} Hz")
    print(f"Total frequency ratio from all intervals: {total_ratio.numerator}/{total_ratio.denominator}")
    print("\nFinal Equation:")
    print(f"{start_freq.numerator}/{start_freq.denominator} Hz * {total_ratio.numerator}/{total_ratio.denominator} = {a}/{b} Hz")
    
    print("\nThe final frequency as a fraction a/b, where a and b are coprime, is given by the list [a, b]:")
    print(f"[{a},{b}]")

calculate_final_frequency()
<<<[62720, 243]>>>