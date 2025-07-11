import math
from fractions import Fraction

def solve_haenschen_klein():
    """
    Solves the problem of finding the final frequency of "Hänschen klein"
    when played with chained just intonation intervals.
    """
    # Define standard just intonation ratios for the intervals found in the melody.
    ratios = {
        "Minor Second up": Fraction(16, 15), "Minor Second down": Fraction(15, 16),
        "Major Second up": Fraction(9, 8),   "Major Second down": Fraction(8, 9),
        "Minor Third up":  Fraction(6, 5),   "Minor Third down":  Fraction(5, 6),
        "Major Third up":  Fraction(5, 4),   "Major Third down":  Fraction(4, 5),
        "Perfect Fifth down": Fraction(2, 3)
    }

    # Define the sequence of musical intervals for the melody "Hänschen klein".
    interval_sequence = [
        ratios["Minor Third down"],  # G -> E
        ratios["Minor Second up"],   # E -> F
        ratios["Minor Third down"],  # F -> D
        ratios["Major Second down"], # D -> C
        ratios["Major Second up"],   # C -> D
        ratios["Major Second up"],   # D -> E
        ratios["Minor Second up"],   # E -> F
        ratios["Major Second up"],   # F -> G
        ratios["Minor Third down"],  # G -> E
        ratios["Minor Second up"],   # E -> F
        ratios["Minor Third down"],  # F -> D
        ratios["Major Second down"], # D -> C
        ratios["Major Third up"],    # C -> E
        ratios["Minor Third up"],    # E -> G
        ratios["Perfect Fifth down"],# G -> C
        ratios["Minor Third down"],  # C -> A
        ratios["Major Second down"], # A -> G
        ratios["Major Second down"], # G -> F
        ratios["Minor Second down"], # F -> E
        ratios["Major Second down"], # E -> D
        ratios["Major Second up"],   # D -> E
        ratios["Minor Second up"],   # E -> F
        ratios["Major Second up"],   # F -> G
        ratios["Major Second up"],   # G -> A
        ratios["Major Second down"], # A -> G
        ratios["Major Second down"], # G -> F
        ratios["Minor Second down"], # F -> E
        ratios["Minor Second up"],   # E -> F
        ratios["Major Second up"],   # F -> G
        ratios["Major Second up"],   # G -> A
        ratios["Major Second up"],   # A -> B
        ratios["Minor Second up"],   # B -> C'
    ]

    # Start with the initial frequency of G = 392 Hz
    final_freq = Fraction(392, 1)
    
    # Build the full equation string while calculating the final frequency
    equation_string = f"{final_freq.numerator}/{final_freq.denominator}"
    
    for r in interval_sequence:
        final_freq *= r
        equation_string += f" * {r.numerator}/{r.denominator}"

    equation_string += f" = {final_freq.numerator}/{final_freq.denominator}"

    # Print the full calculation as a single equation
    print("The final frequency is calculated by chained multiplication of just interval ratios.")
    print("Full calculation:")
    print(equation_string)

    # The result is the numerator and denominator of the final frequency
    a = final_freq.numerator
    b = final_freq.denominator

    print("\nFinal frequency in a/b Hertz form [a, b]:")
    print([a, b])

solve_haenschen_klein()