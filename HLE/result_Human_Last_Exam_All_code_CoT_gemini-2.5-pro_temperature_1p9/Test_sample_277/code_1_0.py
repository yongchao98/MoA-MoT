import math
from fractions import Fraction

def calculate_final_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein"
    using path-dependent just intonation.
    """
    
    # Step 1: Define the starting frequency of the first note (G).
    start_freq = Fraction(392, 1)

    # Step 2: Define the sequence of just interval ratios for the entire verse.
    # The melody is: G E E F D D C D E F G G G | G E E F D D C E G G C C C
    # The transitions and their corresponding ratios are listed below.
    interval_ratios = [
        # --- First Line ---
        Fraction(4, 5),   # 1. G -> E (Descending Major Third)
        Fraction(1, 1),   # 2. E -> E (Unison)
        Fraction(16, 15), # 3. E -> F (Ascending Minor Second)
        Fraction(5, 6),   # 4. F -> D (Descending Minor Third)
        Fraction(1, 1),   # 5. D -> D (Unison)
        Fraction(8, 9),   # 6. D -> C (Descending Major Second)
        Fraction(9, 8),   # 7. C -> D (Ascending Major Second)
        Fraction(9, 8),   # 8. D -> E (Ascending Major Second)
        Fraction(16, 15), # 9. E -> F (Ascending Minor Second)
        Fraction(9, 8),   # 10. F -> G (Ascending Major Second)
        Fraction(1, 1),   # 11. G -> G (Unison)
        Fraction(1, 1),   # 12. G -> G (Unison)
        # --- Transition between lines ---
        Fraction(1, 1),   # 13. G -> G (Unison)
        # --- Second Line ---
        Fraction(4, 5),   # 14. G -> E (Descending Major Third)
        Fraction(1, 1),   # 15. E -> E (Unison)
        Fraction(16, 15), # 16. E -> F (Ascending Minor Second)
        Fraction(5, 6),   # 17. F -> D (Descending Minor Third)
        Fraction(1, 1),   # 18. D -> D (Unison)
        Fraction(8, 9),   # 19. D -> C (Descending Major Second)
        Fraction(5, 4),   # 20. C -> E (Ascending Major Third)
        Fraction(6, 5),   # 21. E -> G (Ascending Minor Third)
        Fraction(1, 1),   # 22. G -> G (Unison)
        Fraction(2, 3),   # 23. G -> C (Descending Perfect Fifth)
        Fraction(1, 1),   # 24. C -> C (Unison)
        Fraction(1, 1),   # 25. C -> C (Unison)
    ]

    # Step 3: Calculate the total ratio by multiplying all individual ratios.
    total_ratio = Fraction(1, 1)
    for ratio in interval_ratios:
        total_ratio *= ratio

    # Step 4: Calculate the final frequency by applying the total ratio to the starting frequency.
    final_freq = start_freq * total_ratio
    
    # Step 5: Extract the coprime numerator (a) and denominator (b).
    a = final_freq.numerator
    b = final_freq.denominator

    # The problem requires outputting the numbers in the final equation.
    # The final equation is: start_freq * total_ratio = final_freq
    # Which can be written as: (392/1) * (total_ratio_num/total_ratio_den) = a/b
    print(f"Starting Frequency: {start_freq.numerator}/{start_freq.denominator} Hz")
    print(f"Total Ratio Multiplier: {total_ratio.numerator}/{total_ratio.denominator}")
    print(f"Final Frequency Equation: ({start_freq.numerator}/{start_freq.denominator}) * ({total_ratio.numerator}/{total_ratio.denominator}) = {a}/{b} Hz")
    
    # Finally, print the result in the requested list format.
    print(f"The resulting list [a, b] is:")
    print([a, b])

calculate_final_frequency()