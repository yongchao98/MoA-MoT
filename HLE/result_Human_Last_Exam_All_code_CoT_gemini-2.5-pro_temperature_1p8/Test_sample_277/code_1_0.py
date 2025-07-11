from fractions import Fraction

def calculate_final_frequency():
    """
    Calculates the frequency of the last note of "Hänschen klein" using just intonation.

    The function follows these steps:
    1. Sets the starting frequency and the note sequence for a full verse.
    2. Defines the frequency ratios for intervals in a just-intonation G-major scale.
    3. Iterates through the melody, calculating the cumulative frequency ratio.
    4. Computes the final frequency and presents the result as requested.
    """
    # 1. Setup initial frequency and the melody for one full verse
    # ending with "...läuft nach Haus geschwind."
    start_freq = Fraction(392, 1)

    # The melody combines several standard phrases to cover the full verse text.
    # Phrase 1&2: "Hänschen klein ... wohlgemut." (G A B C D D, D C B A G G) x2
    # Phrase 3: "Aber Mutter weinet sehr ..." (D D E F# G G, G F# E D C C)
    # Phrase 4: "Da besinnt sich das Kind ..." (C B A G A B, C C B A G G)
    notes_str = "G A B C D D D C B A G G G A B C D D D C B A G G D D E Fsharp G G G Fsharp E D C C C B A G A B C C B A G G"
    melody = notes_str.replace("Fsharp", "F#").split()

    # 2. Define just intonation ratios for G-major (5-limit tuning)
    # The scale steps are: Major Tone (T=9/8), Minor Tone (t=10/9), Semitone (S=16/15)
    # G-major pattern: G --T-- A --t-- B --S-- C --T-- D --t-- E --T-- F# --S-- G
    ratios = {
        ('G', 'A'): Fraction(9, 8),
        ('A', 'B'): Fraction(10, 9),
        ('B', 'C'): Fraction(16, 15),
        ('C', 'D'): Fraction(9, 8),
        ('D', 'E'): Fraction(10, 9),
        ('E', 'F#'): Fraction(9, 8),
        ('F#', 'G'): Fraction(16, 15),
    }
    # Add ratios for descending steps
    for (note1, note2), ratio in list(ratios.items()):
        ratios[(note2, note1)] = 1 / ratio

    # Add the special ratio for the melodic jump G -> D (Perfect Fifth)
    ratios[('G', 'D')] = Fraction(3, 2)
    
    # 3. Calculate the cumulative frequency ratio by iterating through the melody
    total_ratio = Fraction(1, 1)
    for i in range(len(melody) - 1):
        note1 = melody[i]
        note2 = melody[i+1]

        if note1 == note2:
            continue  # Ratio is 1, no change
        
        interval_ratio = ratios.get((note1, note2))
        if interval_ratio is None:
            raise ValueError(f"The interval ratio for {note1} -> {note2} is not defined.")
            
        total_ratio *= interval_ratio

    # 4. Calculate the final frequency and extract the numerator and denominator
    final_freq = start_freq * total_ratio
    a = final_freq.numerator
    b = final_freq.denominator
    
    # The problem requests to output the numbers in the final equation.
    # The final equation is: F_final = F_initial * Total_Ratio
    print(f"The calculation is based on the following equation:")
    print(f"Final Frequency = Initial Frequency * Total Ratio")
    print(f"({a}/{b}) Hz = ({start_freq.numerator}/{start_freq.denominator}) Hz * ({total_ratio.numerator}/{total_ratio.denominator})")
    
    # Final answer format
    print(f"\nThe frequency of the last played note is {a}/{b} Hertz.")
    print("The requested list [a, b] is:")
    print([a, b])

calculate_final_frequency()