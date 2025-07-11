import fractions

def solve_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    using just intonation intervals chained together.
    """

    # Define the simple integer ratios for just (pure) intervals.
    # The ratio for a descending interval is the inverse of the ascending one.
    ratios = {
        'm2_up': fractions.Fraction(16, 15), 'm2_down': fractions.Fraction(15, 16),
        'M2_up': fractions.Fraction(9, 8),   'M2_down': fractions.Fraction(8, 9),
        'm3_up': fractions.Fraction(6, 5),   'm3_down': fractions.Fraction(5, 6),
        'M3_up': fractions.Fraction(5, 4),   'M3_down': fractions.Fraction(4, 5),
        'P4_up': fractions.Fraction(4, 3),   'P4_down': fractions.Fraction(3, 4),
        'P5_up': fractions.Fraction(3, 2),   'P5_down': fractions.Fraction(2, 3),
    }

    # This is the sequence of musical intervals for one verse of the song.
    # Unison intervals (note to same note) have a ratio of 1 and are omitted.
    # Melody: G E E | F D D | C D E F | G G G | G E E | F D D | C E G G | C C C
    interval_names = [
        'm3_down',  # G -> E
        'm2_up',    # E -> F
        'M3_down',  # F -> D
        'M2_down',  # D -> C
        'M2_up',    # C -> D
        'M2_up',    # D -> E
        'm2_up',    # E -> F
        'M2_up',    # F -> G
        # First musical phrase ends on G. Second phrase starts.
        'm3_down',  # G -> E
        'm2_up',    # E -> F
        'M3_down',  # F -> D
        'M2_down',  # D -> C
        'M3_up',    # C -> E
        'm3_up',    # E -> G
        'P5_down',  # G -> C
        # The song ends on C, held for 3 beats.
    ]

    # The starting frequency of the first note (G) is 392 Hz.
    f_start = fractions.Fraction(392)
    
    # Calculate the total ratio by multiplying the ratios of all intervals.
    total_ratio = fractions.Fraction(1, 1)
    for name in interval_names:
        total_ratio *= ratios[name]
        
    # The final frequency is the starting frequency multiplied by the total ratio.
    f_final = f_start * total_ratio
    
    # The fractions module automatically simplifies the fraction to its coprime numerator and denominator.
    a = f_final.numerator
    b = f_final.denominator
    
    # Print the equation and the final answer in the requested format.
    print(f"Starting Frequency: {f_start.numerator}/{f_start.denominator} Hz")
    print(f"Total Ratio: {total_ratio.numerator}/{total_ratio.denominator}")
    print(f"Final Frequency = ({f_start.numerator} * {total_ratio.numerator}) / {total_ratio.denominator} = {a}/{b} Hz")
    print("The final answer is [a, b]:")
    print(f"[{a},{b}]")

solve_frequency()