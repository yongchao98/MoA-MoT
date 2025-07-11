import sys

def solve_key_signature_formula():
    """
    This function derives the formula for the sum of sharps as requested.
    It calculates the base sum (S_0) and the coefficient for n,
    then prints the derivation and the final formula.
    """
    
    # Step 1: Calculate the base sum for n=0 (S_0).
    # The 12 notes are C, C#, D, D#, E, F, F#, G, G#, A, A#, B.
    # We list the number of sharps in each major key.
    # For keys with flats, we convert using the rule: sharps = 12 - flats.
    
    # Note: Sharps/Flats for standard major keys
    # C: 0 sharps
    # C#: 7 sharps
    # D: 2 sharps
    # D# (enharmonic to Eb): 3 flats -> 12 - 3 = 9 sharps
    # E: 4 sharps
    # F: 1 flat -> 12 - 1 = 11 sharps
    # F#: 6 sharps
    # G: 1 sharp
    # G# (enharmonic to Ab): 4 flats -> 12 - 4 = 8 sharps
    # A: 3 sharps
    # A# (enharmonic to Bb): 2 flats -> 12 - 2 = 10 sharps
    # B: 5 sharps
    
    base_sharps = {
        'C': 0, 'C#': 7, 'D': 2, 'D#': 9, 'E': 4, 'F': 11,
        'F#': 6, 'G': 1, 'G#': 8, 'A': 3, 'A#': 10, 'B': 5
    }
    
    s0 = sum(base_sharps.values())
    
    # Step 2: Determine the effect of sharping the tonics n times.
    # Adding one sharp to a key's tonic adds 7 sharps to its key signature.
    # This must be done for all 12 notes.
    
    sharps_per_n_increment = 12 * 7
    
    # Step 3: Construct and print the final formula S_n = S_0 + (12 * 7)n
    
    print("Derivation of the formula S_n for the sum of sharps:")
    print("-" * 50)
    
    print("1. Calculate the sum of sharps for the base case (n=0), S_0.")
    print("   The number of sharps for each of the 12 major keys (with flats converted to sharps) are:")
    print(f"   {list(base_sharps.values())}")
    print(f"   The sum, S_0 = { ' + '.join(map(str, base_sharps.values())) } = {s0}")
    
    print("\n2. Calculate the number of sharps added for a given 'n'.")
    print("   Adding one sharp to a tonic adds 7 sharps to its key signature.")
    print("   Since this is done for all 12 notes 'n' times, the total sharps added is:")
    print(f"   12 notes * 7 sharps/note * n = {sharps_per_n_increment}n")
    
    print("\n3. Combine the terms to find the final formula: S_n = S_0 + 84n")
    print("   Substituting the calculated value of S_0, the final equation is:")
    
    # Final formatted output of the equation
    # The user requested to "output each number in the final equation"
    sys.stdout.write("S_n = ")
    sys.stdout.flush() # Flush to ensure print order
    # Print the numbers individually
    print(s0, end=' + ')
    print(sharps_per_n_increment, end='n\n')

solve_key_signature_formula()