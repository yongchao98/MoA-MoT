from fractions import Fraction

def solve_hanschen_klein_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation.
    """
    
    # The starting frequency of the first note (G) is 392 Hz.
    # We represent it as a Fraction for precise arithmetic.
    start_frequency = Fraction(392, 1)

    # Define the just intonation ratios for the musical intervals present in the song.
    # In just intonation, the size of a whole tone can differ based on its position in the scale.
    # D->E is a "minor tone" (10/9) in a C-major context, while others are "major tones" (9/8).
    ratios_map = {
        "unison": Fraction(1, 1),
        "min2_asc": Fraction(16, 15), # E.g., E -> F
        "maj2_asc_major": Fraction(9, 8), # "Major tone", e.g., C -> D, F -> G
        "maj2_asc_minor": Fraction(10, 9), # "Minor tone", e.g., D -> E
        "maj2_desc": Fraction(8, 9), # E.g., D -> C
        "min3_asc": Fraction(6, 5), # E.g., E -> G
        "min3_desc": Fraction(5, 6), # E.g., F -> D
        "maj3_asc": Fraction(5, 4), # E.g., C -> E
        "maj3_desc": Fraction(4, 5), # E.g., G -> E
        "p5_desc": Fraction(2, 3), # E.g., G -> C
    }

    # The melody consists of 26 notes, which means there are 25 transitions.
    # We list the type of interval for each transition.
    # Melody pieces: A(GEEFDD), B(CDEFGGG), A'(GEEFDD), C(CEGGCCC)
    transitions = [
        # A: G -> E -> E -> F -> D -> D
        "maj3_desc", "unison", "min2_asc", "min3_desc", "unison",
        # Transition A -> B: D -> C
        "maj2_desc",
        # B: C -> D -> E -> F -> G -> G -> G
        "maj2_asc_major", "maj2_asc_minor", "min2_asc", "maj2_asc_major", "unison", "unison",
        # Transition B -> A': G -> G
        "unison",
        # A': G -> E -> E -> F -> D -> D
        "maj3_desc", "unison", "min2_asc", "min3_desc", "unison",
        # Transition A' -> C: D -> C
        "maj2_desc",
        # C: C -> E -> G -> G -> C -> C -> C
        "maj3_asc", "min3_asc", "unison", "p5_desc", "unison", "unison",
    ]

    # Calculate the cumulative ratio by multiplying the ratios of all transitions.
    total_ratio = Fraction(1, 1)
    for t_name in transitions:
        total_ratio *= ratios_map[t_name]
        
    # The final frequency is the starting frequency multiplied by the total ratio.
    final_frequency = start_frequency * total_ratio

    # The result is a/b Hertz, where a and b are coprime.
    a = final_frequency.numerator
    b = final_frequency.denominator

    print(f"The starting frequency of G is: {start_frequency.numerator}/{start_frequency.denominator} Hz")
    print(f"The cumulative frequency ratio after one verse is: {total_ratio.numerator}/{total_ratio.denominator}")
    print(f"The frequency of the last note is the product: ({start_frequency.numerator}/{start_frequency.denominator}) * ({total_ratio.numerator}/{total_ratio.denominator})")
    print(f"This results in a final frequency of: {a}/{b} Hertz")
    print("The answer in the format [a, b] is:")
    print([a, b])


solve_hanschen_klein_frequency()