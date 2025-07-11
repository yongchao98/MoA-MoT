import math
from fractions import Fraction

def solve():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    using just intonation intervals.
    """
    # 1. Initial setup
    # The song starts on G tuned to 392 Hz.
    initial_freq = Fraction(392, 1)

    # Standard 5-limit just intonation ratios for relevant intervals
    ratios = {
        "M3_down": Fraction(4, 5),   # Major Third down (e.g., G -> E)
        "m2_up":   Fraction(16, 15), # Minor Second up (e.g., E -> F)
        "M2_up":   Fraction(9, 8),   # Major Second up (e.g., D -> E)
        "M2_down": Fraction(8, 9),   # Major Second down (e.g., E -> D)
        "P5_up":   Fraction(3, 2),   # Perfect Fifth up (e.g., C -> G)
        "m2_down": Fraction(15, 16), # Minor Second down (e.g., F -> E)
        "P4_up":   Fraction(4, 3),   # Perfect Fourth up (e.g., D -> G)
        "P5_down": Fraction(2, 3),   # Perfect Fifth down (e.g., G -> C)
    }

    # 2. Calculate the ratio for melodic line A: G E E F D D E D C G G G
    # The intervals are: G->E, E->F, F->D, D->E, E->D, D->C, C->G
    # Unison intervals (E->E, D->D, G->G) have a ratio of 1.
    line_A_intervals = [
        ratios["M3_down"],
        ratios["m2_up"],
        ratios["M3_down"],
        ratios["M2_up"],
        ratios["M2_down"],
        ratios["M2_down"],
        ratios["P5_up"],
    ]
    # The overall ratio for line A is the product of these interval ratios.
    R_A = math.prod(line_A_intervals)

    # 3. Calculate the ratio for playing melodic line B after line A.
    # Last note of A is G. First note of B is A. Melody: A A G F F E D G C C C C
    # The intervals are: G->A, A->G, G->F, F->E, E->D, D->G, G->C
    line_B_full_intervals = [
        ratios["M2_up"],    # Transition G->A
        ratios["M2_down"],
        ratios["M2_down"],
        ratios["m2_down"],
        ratios["M2_down"],
        ratios["P4_up"],
        ratios["P5_down"],
    ]
    R_B = math.prod(line_B_full_intervals)
    
    # 4. Calculate the ratio for the transition from line B back to line A.
    # Last note of B is C. First note of A is G.
    trans_B_to_A = ratios["P5_up"]

    # 5. Assemble the final calculation based on the AABA structure.
    # The overall ratio is a result of playing line A, then A again,
    # then B, then transitioning to A and playing it a final time.
    # Note: The G->G transition between the first two A-lines is a unison (ratio 1).
    total_ratio = (R_A ** 3) * R_B * trans_B_to_A
    
    final_freq = initial_freq * total_ratio

    # 6. Output the results
    print("This script calculates the final frequency of 'Hänschen klein'.")
    print("\nMelody structure: AABA")
    print(f"Starting frequency of G: {initial_freq.numerator}/{initial_freq.denominator} Hz")
    
    print("\nComponent Ratios:")
    print(f"Ratio for line A (G E E...): {R_A.numerator}/{R_A.denominator}")
    print(f"Ratio for playing line B after A (G -> A A G...): {R_B.numerator}/{R_B.denominator}")
    print(f"Ratio for transition from line B to A (C->G): {trans_B_to_A.numerator}/{trans_B_to_A.denominator}")

    print("\nFinal Equation:")
    print(f"Final Frequency = Initial_Freq * (R_A)^3 * R_B * R_trans_B_to_A")
    print(f"Final Frequency = ({initial_freq}) * ({R_A.numerator}/{R_A.denominator})^3 * ({R_B.numerator}/{R_B.denominator}) * ({trans_B_to_A.numerator}/{trans_B_to_A.denominator})")

    a = final_freq.numerator
    b = final_freq.denominator
    
    print("\nEach number in the final equation:")
    print(f"Initial Frequency: {initial_freq.numerator}")
    print(f"Line A Ratio: {R_A.numerator}, {R_A.denominator} (used 3 times)")
    print(f"Line B Full Ratio: {R_B.numerator}, {R_B.denominator} (used 1 time)")
    print(f"Transition B->A Ratio: {trans_B_to_A.numerator}, {trans_B_to_A.denominator} (used 1 time)")
    
    print(f"\nFinal frequency = {a} / {b} Hertz.")
    print("\nThe result [a, b] is:")
    print([a, b])
    
    # Return the final answer in the specified format for the platform
    return [a, b]

result = solve()
# The final result is extracted and formatted below by the system.
# The `solve()` function prints the detailed explanation.
# For example: <<< [6734508720128, 2562890625] >>>
final_answer_str = f"<<<[{result[0]}, {result[1]}]>>>"
# This print is for display purposes, the actual answer is captured from the return value.
# print(final_answer_str)