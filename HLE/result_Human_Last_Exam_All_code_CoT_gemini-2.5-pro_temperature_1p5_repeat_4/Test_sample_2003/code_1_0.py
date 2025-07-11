def solve_key_signature_formula():
    """
    Derives and prints a formula for the sum of sharps in key signatures.
    """
    # 1. Represent the 12 initial notes as a list of tuples: (NoteLetter, initial_sharps)
    initial_notes = [
        ('C', 0), ('C', 1),  # C, C#
        ('D', 0), ('D', 1),  # D, D#
        ('E', 0),            # E
        ('F', 0), ('F', 1),  # F, F#
        ('G', 0), ('G', 1),  # G, G#
        ('A', 0), ('A', 1),  # A, A#
        ('B', 0)             # B
    ]

    # 2. Map natural note letters to the number of sharps in their standard key signature.
    # Flats are represented as negative sharps (e.g., F major has 1 flat = -1 sharps).
    base_key_sharps = {
        'C': 0, 'D': 2, 'E': 4, 'F': -1, 'G': 1, 'A': 3, 'B': 5
    }

    # 3. For the formula S(n) = C + D*n, calculate the components.
    # The total sum is Sum[base_sharps(note_i) + 7 * (s_i + n)] over all 12 notes.
    # This can be expanded to: Sum[base_sharps] + 7*Sum[s_i] + 12*7*n

    # Calculate Sum[base_sharps(note_i)]
    sum_base_sharps = 0
    for note_letter, _ in initial_notes:
        sum_base_sharps += base_key_sharps[note_letter]

    # Calculate Sum[s_i]
    sum_initial_sharps = 0
    for _, s_i in initial_notes:
        sum_initial_sharps += s_i
        
    # 4. As analyzed, for n > 0, the number of sharps is never negative,
    # so we don't need to apply the flat-to-sharp conversion rule.
    # The rule for F major -> E sharp major (11 sharps) only matters for the n=0 case.

    # 5. Calculate the constant part 'C' of the formula S(n) = C + D*n
    # C = Sum[base_sharps] + 7*Sum[s_i]
    constant_term = sum_base_sharps + 7 * sum_initial_sharps

    # The 'n' coefficient 'D' is 12 notes * 7 sharps per full cycle.
    n_coefficient = 12 * 7

    # Print the derived formula and the numbers used in it.
    print("The formula for the sum of sharps S(n) for n > 0 is derived as follows:")
    print(f"Sum = (Sum of base key sharps) + 7 * (Sum of initial sharps on notes) + n * (12 notes * 7)")
    print(f"Sum = ({sum_base_sharps}) + 7 * ({sum_initial_sharps}) + n * ({n_coefficient})")
    print(f"Sum = ({sum_base_sharps}) + ({7 * sum_initial_sharps}) + {n_coefficient}n")
    print(f"Sum = {constant_term} + {n_coefficient}n")
    
    print("\nSimplified Formula:")
    print(f"S(n) = {constant_term} + {n_coefficient}n")

solve_key_signature_formula()
<<<54 + 84n>>>