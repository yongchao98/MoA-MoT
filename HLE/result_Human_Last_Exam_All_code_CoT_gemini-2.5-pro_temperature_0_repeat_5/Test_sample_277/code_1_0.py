from fractions import Fraction

def calculate_final_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation.

    The calculation is based on a full verse of the song, starting with G=392 Hz.
    The frequency of each note is determined by multiplying the frequency of the
    preceding note by the just intonation ratio of the interval between them.
    """

    # Initial frequency of the first note (G)
    initial_freq = Fraction(392, 1)

    # The melody for a full verse is composed of two main phrases.
    # Phrase 1: G E E F D D C D E F G G G
    # Phrase 2: G E E F D D C E G G C C C
    # The last note is C.

    # Ratios for the intervals in the first musical phrase.
    # Intervals: G->E, E->E, E->F, F->D, D->D, D->C, C->D, D->E, E->F, F->G, G->G, G->G
    ratios_phrase1 = [
        Fraction(5, 6),   # Minor Third down
        Fraction(1, 1),   # Unison
        Fraction(16, 15), # Minor Second up
        Fraction(5, 6),   # Minor Third down
        Fraction(1, 1),   # Unison
        Fraction(8, 9),   # Major Second down
        Fraction(9, 8),   # Major Second up
        Fraction(9, 8),   # Major Second up
        Fraction(16, 15), # Minor Second up
        Fraction(9, 8),   # Major Second up
        Fraction(1, 1),   # Unison
        Fraction(1, 1),   # Unison
    ]

    # Ratios for the intervals in the second musical phrase, starting from the last note of phrase 1.
    # Intervals: G->G, G->E, E->E, E->F, F->D, D->D, D->C, C->E, E->G, G->G, G->C, C->C, C->C
    ratios_phrase2 = [
        Fraction(1, 1),   # Unison
        Fraction(5, 6),   # Minor Third down
        Fraction(1, 1),   # Unison
        Fraction(16, 15), # Minor Second up
        Fraction(5, 6),   # Minor Third down
        Fraction(1, 1),   # Unison
        Fraction(8, 9),   # Major Second down
        Fraction(5, 4),   # Major Third up
        Fraction(6, 5),   # Minor Third up
        Fraction(1, 1),   # Unison
        Fraction(2, 3),   # Perfect Fifth down
        Fraction(1, 1),   # Unison
        Fraction(1, 1),   # Unison
    ]

    all_ratios = ratios_phrase1 + ratios_phrase2

    # Calculate the final frequency by multiplying all ratios
    current_freq = initial_freq
    equation_parts = [str(initial_freq.numerator)]

    for ratio in all_ratios:
        current_freq *= ratio
        # Add the ratio to the equation string, skipping unisons for brevity
        if ratio != 1:
            equation_parts.append(f"{ratio.numerator}/{ratio.denominator}")

    final_freq = current_freq
    a = final_freq.numerator
    b = final_freq.denominator

    # Print the full equation
    print("The calculation for the final frequency is:")
    equation_str = " * ".join(equation_parts)
    print(f"{equation_str} = {a}/{b}")
    
    # Print the final answer in the required format
    print("\nThe frequency of the last played note is a/b Hertz, where [a, b] is:")
    print(f"[{a}, {b}]")
    
    return [a, b]

# Execute the calculation and store the final answer
final_answer = calculate_final_frequency()

# The final answer in the specified format
# print(f"\n<<<{final_answer}>>>")