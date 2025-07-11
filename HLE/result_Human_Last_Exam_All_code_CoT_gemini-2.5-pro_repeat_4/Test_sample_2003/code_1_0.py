def solve_key_signature_sum():
    """
    This script calculates the sum of sharps across all 12 major key signatures
    after a global transformation of sharping every note 'n' times.
    """
    # We choose an arbitrary n > 0 for demonstration, as per the problem.
    # The logic will show the result is independent of n.
    n = 1
    
    print(f"This analysis demonstrates the calculation for n = {n}.")
    print("The final formula, however, will be constant for any n > 0.\n")

    # Step 1: Represent the 12 chromatic pitches as integers 0-11.
    initial_pitches = list(range(12))

    # Step 2: Apply the transformation. "Sharping n times" means adding n to the pitch value.
    # The operation is cyclical (modulo 12).
    # The resulting set of pitches is just a permutation of the original set.
    transformed_pitches = [(p + n) % 12 for p in initial_pitches]

    # Step 3: Calculate the number of sharps for each new key signature.
    # The number of sharps for a major key with a tonic of pitch 'p' is (p * 7) % 12.
    # This formula is derived from the circle of fifths and correctly converts
    # flat keys to their enharmonic sharp equivalents as requested.
    sharps_in_each_key = []
    for p in transformed_pitches:
        num_sharps = (p * 7) % 12
        sharps_in_each_key.append(num_sharps)

    # Step 4: Sum the sharps.
    # Since `transformed_pitches` is a permutation of {0, ..., 11} and multiplying by 7 (which is
    # coprime to 12) also creates a permutation, the list `sharps_in_each_key` will always
    # contain the numbers {0, 1, ..., 11} in some order.
    # We sort it here just for a clean, ordered display of the equation.
    sharps_in_each_key.sort()
    
    # Calculate the total sum.
    total_sharps = sum(sharps_in_each_key)

    # Step 5: Display the final equation and the derived formula.
    # As requested, we output each number in the final sum.
    equation_str = " + ".join(map(str, sharps_in_each_key))
    print("The number of sharps for the 12 key signatures will always be a permutation of the set {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}.")
    print("\nThe final sum is calculated as follows:")
    print(f"{equation_str} = {total_sharps}")
    
    print("\n---")
    print("The sum is constant because the set of notes remains the same after the transformation.")
    print("The sum is always the sum of integers from 0 to 11, which is 66.")
    print("Therefore, the simplified formula for the sum S(n) is a constant.")
    print("S(n) = 66")

solve_key_signature_sum()