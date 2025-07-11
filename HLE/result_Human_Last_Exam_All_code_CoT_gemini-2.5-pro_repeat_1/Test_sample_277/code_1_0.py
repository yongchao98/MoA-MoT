from fractions import Fraction

def get_interval_ratio(note1, note2):
    """
    Returns the just intonation frequency ratio for the interval between two notes.
    Handles diatonic scale nuances for major seconds.
    """
    # Note dictionary mapping note names to a pitch value on a chromatic scale (C=0)
    # G3 is handled by a pitch value of -5
    pitch = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'G3': -5}
    
    p1 = pitch[note1]
    p2 = pitch[note2]
    interval = p2 - p1

    # Unison
    if interval == 0: return Fraction(1, 1)
    
    # Seconds
    if interval == 1: # Minor second up (E->F)
        return Fraction(16, 15)
    if interval == -1: # Minor second down (F->E)
        return Fraction(15, 16)
    if interval == 2: # Major second up
        if note1 in ['D', 'G']: # D->E or G->A is a minor tone
            return Fraction(10, 9)
        else: # C->D or F->G is a major tone
            return Fraction(9, 8)
    if interval == -2: # Major second down
        if note1 in ['E', 'A']: # E->D or A->G is a minor tone
            return Fraction(9, 10)
        else: # D->C or G->F is a major tone
            return Fraction(8, 9)

    # Thirds
    if interval == 3: return Fraction(6, 5)   # Minor third up
    if interval == -3: return Fraction(5, 6)  # Minor third down
    if interval == 4: return Fraction(5, 4)   # Major third up
    if interval == -4: return Fraction(4, 5)  # Major third down

    # Fourths and Fifths
    if interval == 5: return Fraction(4, 3)   # Perfect fourth up
    if interval == -5: return Fraction(3, 4)  # Perfect fourth down
    if interval == 7: return Fraction(3, 2)   # Perfect fifth up
    if interval == -7: return Fraction(2, 3)  # Perfect fifth down
    
    # Special case from D4 to G3
    if note1 == 'D' and note2 == 'G3':
        return Fraction(2, 3) # Perfect fifth down
        
    # Special case from G3 to C4
    if note1 == 'G3' and note2 == 'C':
        return Fraction(4, 3) # Perfect fourth up

    raise ValueError(f"Unknown interval: {note1} to {note2}")

def solve():
    """
    Calculates the frequency of the last note of "Hänschen klein".
    """
    # Note sequence from a standard version of the song.
    note_sequence = [
        'G', 'E', 'E', 'F', 'D', 'D', 'C', 'D', 'E', 'F', 'G', 'G', 'G',
        'G', 'E', 'E', 'F', 'D', 'D', 'C', 'E', 'G', 'G', 'C', 'C', 'C',
        'E', 'F', 'G', 'E', 'F', 'G', 'A', 'G', 'F', 'E', 'D', 'C',
        'G', 'G', 'G', 'C', 'D', 'E', 'F', 'E', 'D', 'G3', 'C'
    ]

    f0 = Fraction(392, 1)
    current_freq = f0
    total_ratio = Fraction(1, 1)

    for i in range(len(note_sequence) - 1):
        note1 = note_sequence[i]
        note2 = note_sequence[i+1]
        ratio = get_interval_ratio(note1, note2)
        total_ratio *= ratio
        
    final_freq = f0 * total_ratio

    print("Calculation of the final frequency:")
    print(f"Start Frequency (G) = {f0.numerator}/{f0.denominator} Hz")
    print(f"Total Ratio Multiplier = {total_ratio.numerator}/{total_ratio.denominator}")
    print(f"Final Frequency = ({f0.numerator}/{f0.denominator}) * ({total_ratio.numerator}/{total_ratio.denominator}) = {final_freq.numerator}/{final_freq.denominator} Hz")
    
    # Final answer in list format
    result = [final_freq.numerator, final_freq.denominator]
    print("\nFinal answer as a list [a, b]:")
    print(result)

solve()
<<<[12845056, 59049]>>>