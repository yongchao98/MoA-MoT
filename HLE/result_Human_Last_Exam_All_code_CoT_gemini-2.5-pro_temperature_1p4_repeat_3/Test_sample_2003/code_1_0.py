def derive_sharps_formula(n):
    """
    This function calculates the sum of sharps for the 12 major keys
    after each of the 12 chromatic notes has been sharpened n times.

    Args:
        n (int): The number of times to sharp each note (must be > 0).
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # Step 1: Define the initial 12 pitch classes.
    initial_pitch_classes = list(range(12))

    # Step 2: Apply the 'n' transformation.
    # The new set of tonics will have these pitch classes.
    # Note that this is always a permutation of {0, 1, ..., 11}.
    new_tonics_pitch_classes = [(p + n) % 12 for p in initial_pitch_classes]

    # Step 3: Calculate the number of sharps for each new tonic.
    # The formula is (7 * p) mod 12, where p is the pitch class.
    sharps_per_key = [(7 * p) % 12 for p in new_tonics_pitch_classes]

    # Step 4: Calculate the sum.
    total_sharps = sum(sharps_per_key)

    # Step 5: Explain and display the formula and result.
    print(f"For n = {n}, the 12 keys have the following number of sharps:")
    # Sort the list for a clean, predictable output. The sum remains the same.
    sharps_per_key.sort()
    equation_parts = [str(s) for s in sharps_per_key]
    equation_str = " + ".join(equation_parts)
    print(f"Sum = {equation_str} = {total_sharps}")

    print("\nAnalysis:")
    print("The sum is always 66, regardless of the value of n (for n>=0).")
    print("This is because the set of 12 tonic pitch classes is always a permutation of {0, 1, ..., 11},")
    print("and the set of corresponding sharp counts is also always a permutation of {0, 1, ..., 11}.")
    
    # The final simplified formula as an equation with all its numbers.
    final_equation_numbers = sorted([ (7 * p) % 12 for p in range(12) ])
    final_equation_str = " + ".join(map(str, final_equation_numbers))
    
    print("\nThe final simplified formula, S(n), is constant:")
    print(f"S(n) = {final_equation_str} = {sum(final_equation_numbers)}")


# Let's run the derivation for n=1 as an example.
derive_sharps_formula(1)