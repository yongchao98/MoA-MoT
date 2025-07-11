def solve_key_signature_sum():
    """
    Calculates the formula for the sum of sharps in 12 major key signatures
    after each base note has been sharped n times.
    """
    
    # The 12 starting musical notes.
    initial_notes = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    
    # Map note names to a numerical value for calculation (C=0, C#=1, etc.).
    note_to_val = {name: i for i, name in enumerate(initial_notes)}

    # This dictionary defines the "base note" and the number of sharps 's' in the name.
    # Per the problem, F major is treated as E# major, so its base is E and s=1.
    # The base of F# is F.
    note_decomposition = {
        "C":  ("C", 0), "C#": ("C", 1),
        "D":  ("D", 0), "D#": ("D", 1),
        "E":  ("E", 0),
        "F":  ("E", 1),  # F major is rewritten as E# major.
        "F#": ("F", 1),  # F# major's base is F.
        "G":  ("G", 0), "G#": ("G", 1),
        "A":  ("A", 0), "A#": ("A", 1),
        "B":  ("B", 0)
    }

    s0 = 0
    print("Calculating the sum of sharps for n=0 (S₀):")
    print("-" * 45)
    print(f"{'Note':<5} | {'Calculation':<25} | {'Sharps'}")
    print("-" * 45)

    for note in initial_notes:
        base_note_name, s = note_decomposition[note]
        
        # Get the numerical value of the base note.
        base_note_val = note_to_val[base_note_name]
        
        # The number of sharps in the base key signature is (7 * k) % 12.
        # This correctly converts flat keys to their sharp equivalents
        # e.g., F (k=5) -> (7*5)%12 = 11 sharps.
        base_key_sharps = (7 * base_note_val) % 12
        
        # The total sharps for a key is the base key's sharps + 7 for each sharp in the name.
        total_sharps = base_key_sharps + 7 * s
        
        s0 += total_sharps
        
        calc_str = f"({(7 * base_note_val)} % 12) + 7 * {s}"
        print(f"{note:<5} | {calc_str:<25} | {total_sharps}")
        
    print("-" * 45)
    print(f"The total sum for n=0 (S₀) is {s0}.\n")
    
    print("Deriving the formula for S(n):")
    print("Each of the 12 notes is sharped 'n' times.")
    print("This adds 'n' to the sharp count 's' for each key's calculation.")
    print("The increase in sharps for each key is 7 * n.")
    print("The total increase across all 12 keys is 12 * (7 * n) = 84 * n.")
    
    # The coefficient for n is 12 keys * 7 sharps per 'n'
    n_coefficient = 12 * 7
    
    print("\nThe final formula S(n) is S₀ + 84n.")
    print("Substituting the calculated value of S₀:")
    print("\nFinal Equation:")
    print(f"S(n) = {s0} + {n_coefficient}n")

solve_key_signature_sum()
<<<78 + 84n>>>