import math
from fractions import Fraction

# This script calculates the frequency of the last note of "Hänschen klein"
# based on sequential just intonation tuning.

# The melody for a full verse of the Otto Frömmel version is commonly represented as:
# G E E F D D C D E F G G G | G E E F D D C E G G C C C
# This corresponds to the following sequence of intervals between notes.
intervals = [
    'down_minor_third',   # G -> E
    'unison',             # E -> E
    'up_minor_second',    # E -> F
    'down_minor_third',   # F -> D
    'unison',             # D -> D
    'down_major_second',  # D -> C
    'up_major_second',    # C -> D
    'up_major_second',    # D -> E
    'up_minor_second',    # E -> F
    'up_major_second',    # F -> G
    'unison',             # G -> G
    'unison',             # G -> G
    'down_minor_third',   # G -> E
    'unison',             # E -> E
    'up_minor_second',    # E -> F
    'down_minor_third',   # F -> D
    'unison',             # D -> D
    'down_major_second',  # D -> C
    'up_major_third',     # C -> E
    'up_minor_third',     # E -> G
    'unison',             # G -> G
    'down_perfect_fifth', # G -> C
    'unison',             # C -> C
    'unison'              # C -> C (The last note of the verse)
]

# Standard just intonation ratios for musical intervals.
ratios = {
    'unison': Fraction(1, 1),
    'up_minor_second': Fraction(16, 15),
    'down_minor_second': Fraction(15, 16),
    'up_major_second': Fraction(9, 8),
    'down_major_second': Fraction(8, 9),
    'up_minor_third': Fraction(6, 5),
    'down_minor_third': Fraction(5, 6),
    'up_major_third': Fraction(5, 4),
    'down_major_third': Fraction(4, 5),
    'up_perfect_fifth': Fraction(3, 2),
    'down_perfect_fifth': Fraction(2, 3),
}

# The starting frequency is G at 392 Hz. We use fractions for exact arithmetic.
initial_freq = Fraction(392, 1)

# Calculate the total ratio by multiplying all individual interval ratios.
total_ratio = Fraction(1, 1)
for interval in intervals:
    total_ratio *= ratios[interval]

# The final frequency is the initial frequency multiplied by the total ratio.
final_freq = initial_freq * total_ratio

# The result is a fraction a/b, which is already in coprime form
# thanks to the Fraction class.
a = final_freq.numerator
b = final_freq.denominator

# Output the final equation as requested, showing each number.
# Final Frequency = Initial Frequency * (Total Ratio Numerator / Total Ratio Denominator)
print(f"Final equation: {initial_freq.numerator} * {total_ratio.numerator} / {total_ratio.denominator} = {a} / {b} Hertz")

# Output the final answer in the specified list format.
result = [a, b]
print("Result [a, b]:")
print(result)
<<<[25088, 81]>>>