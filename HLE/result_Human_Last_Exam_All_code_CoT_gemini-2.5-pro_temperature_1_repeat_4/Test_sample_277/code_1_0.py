from fractions import Fraction

def solve_song_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation.
    """
    # The plan is to determine the sequence of just intonation interval ratios
    # for the entire song and multiply them together to get a total ratio.
    # This total ratio is then multiplied by the starting frequency.

    # The melody consists of 4 lines. The following list contains the
    # frequency ratios for all 51 intervals in the song, including
    # the transitions between lines. The ratios are derived from a
    # just intonation G-major scale.
    ratios = []

    # Line 1: "Hänschen klein ging allein..." (G E E F# D D C D E F# G G G)
    # Ratios for the 12 intervals in this line, starting from G4.
    ratios.extend([
        Fraction(5, 6), Fraction(1, 1), Fraction(9, 8), Fraction(4, 5), Fraction(1, 1),
        Fraction(8, 9), Fraction(9, 8), Fraction(10, 9), Fraction(9, 8), Fraction(16, 15),
        Fraction(1, 1), Fraction(1, 1)
    ])

    # Line 2: "Stock und Hut stehn ihm gut..." (G E E F# D D C E G G C C C)
    # Transition from L1 end (G4) to L2 start (G4) is a unison.
    ratios.append(Fraction(1, 1))
    ratios.extend([
        Fraction(5, 6), Fraction(1, 1), Fraction(9, 8), Fraction(4, 5), Fraction(1, 1),
        Fraction(8, 9), Fraction(5, 4), Fraction(6, 5), Fraction(1, 1), Fraction(2, 3),
        Fraction(1, 1), Fraction(1, 1)
    ])

    # Line 3: "Aber Mutter weinet sehr..." (A A A B G G F# D D E C C)
    # Transition from L2 end (C4) to L3 start (A4) is a major sixth.
    ratios.append(Fraction(5, 3))
    ratios.extend([
        Fraction(1, 1), Fraction(1, 1), Fraction(10, 9), Fraction(4, 5), Fraction(1, 1),
        Fraction(15, 16), Fraction(4, 5), Fraction(1, 1), Fraction(10, 9), Fraction(4, 5),
        Fraction(1, 1), Fraction(1, 1)
    ])

    # Line 4: "Da besinnt sich das Kind..." (G E E F# D D C E G G C C C)
    # Transition from L3 end (C4) to L4 start (G4) is a perfect fifth.
    ratios.append(Fraction(3, 2))
    ratios.extend([
        Fraction(5, 6), Fraction(1, 1), Fraction(9, 8), Fraction(4, 5), Fraction(1, 1),
        Fraction(8, 9), Fraction(5, 4), Fraction(6, 5), Fraction(1, 1), Fraction(2, 3),
        Fraction(1, 1), Fraction(1, 1)
    ])

    # The starting frequency is G4 = 392 Hz.
    initial_freq = Fraction(392, 1)

    # Calculate the total ratio by multiplying all individual interval ratios.
    total_ratio = Fraction(1, 1)
    for r in ratios:
        total_ratio *= r

    # The final frequency is the initial frequency multiplied by the total ratio.
    final_freq = initial_freq * total_ratio

    # The Fraction object automatically keeps the numerator and denominator coprime.
    a = final_freq.numerator
    b = final_freq.denominator

    print("The frequency of the last note is calculated as follows:")
    print(f"Final Frequency = Initial Frequency * Total Ratio")
    print(f"Final Frequency = {initial_freq.numerator}/{initial_freq.denominator} Hz * {total_ratio.numerator}/{total_ratio.denominator}")
    print(f"Final Frequency = {a}/{b} Hz")
    
    result_list = [a, b]
    print("\nThe answer in the requested list format [a, b] is:")
    print(result_list)

solve_song_frequency()