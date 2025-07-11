def solve_sharps_formula():
    """
    Derives and explains the formula for the sum of sharps S(n).
    The formula is for the sum of the number of sharps of 12 major key signatures
    derived from an initial set of 12 notes, where each note has been sharped n times.
    """
    
    print("Deriving the formula S(n) for the total number of sharps.")
    print("Let n be the number of times the initial 12 notes are sharped (n >= 0).")
    print("-" * 40)

    # Step 1: Calculate the base sum S(0) for n=0
    print("Step 1: Calculate the sum of sharps for the base case (n=0), S(0).")
    print("\nThe initial 12 notes are C, C#, D, D#, E, F, F#, G, G#, A, A#, B.")
    print("Numerically, we can represent these notes' pitch classes as k = 0, 1, ..., 11.")
    print("The number of sharps for a major key whose tonic has pitch class k is given by the formula (7 * k) % 12.")
    print("This formula correctly handles the rule to convert flat keys to sharp keys (e.g., for F major, k=5, (7*5)%12 = 11 sharps for E# major).")

    # Calculate the number of sharps for each of the 12 base keys
    # As k runs from 0 to 11, (7*k)%12 permutes the values 0 to 11 because 7 and 12 are coprime.
    # So the sum is simply the sum of integers from 0 to 11.
    s_0 = sum(range(12))
    
    print(f"\nThe sum for n=0, S(0), is the sum of (7*k)%12 for k=0 to 11, which equals the sum of 0 to 11.")
    print(f"The first number in our final equation is the base sum, S(0) = {s_0}.")
    print("-" * 40)

    # Step 2: Determine the effect of n
    print("Step 2: Analyze the effect of sharping each note n times.")
    print("\nThe rule 'do not simplify key signatures' means adding a sharp to a tonic's name adds 7 sharps to its key signature.")
    print("For example: C major (0 sharps) -> C# major (7 sharps) -> C## major (14 sharps).")
    print("So, sharping a note n times adds (7 * n) sharps to its key signature's count.")
    
    increase_per_note = 7
    total_notes = 12
    coefficient_of_n = increase_per_note * total_notes

    print(f"\nThis increase of (7 * n) sharps applies to all 12 notes.")
    print(f"The total increase in sharps across all notes is 12 * (7 * n).")
    print(f"The second number in our final equation is the coefficient of n, which is 12 * 7 = {coefficient_of_n}.")
    print("-" * 40)

    # Step 3: Combine for the final formula
    print("Step 3: Combine the base sum and the increase to find the final formula S(n).")
    print("\nS(n) = S(0) + (Total increase from n)")
    print("S(n) = (Base Sum) + (Coefficient of n) * n")
    
    print("\nThe final simplified formula is:")
    print(f"Sum = {s_0} + {coefficient_of_n}n")

# Execute the function to print the derivation
solve_sharps_formula()