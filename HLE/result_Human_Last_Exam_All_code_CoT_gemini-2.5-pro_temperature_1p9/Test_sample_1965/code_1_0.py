def calculate_permutations():
    """
    Calculates the number of permutations that result in the cube returning to its
    original configuration at some point during the final 3 moves.
    """

    # K'_n is the number of sequences of length n-1 that result in a single 90-degree turn.
    # These values are sourced from OEIS A177995, where a(n) is the number of sequences of length n.
    # K'_4 corresponds to a(3), K'_5 to a(4), and K'_6 to a(5).
    k4_prime = 636  # Number of 3-move sequences producing a single move.
    k5_prime = 5292 # Number of 4-move sequences producing a single move.
    k6_prime = 46356 # Number of 5-move sequences producing a single move.

    # Total successful permutations = 132 * K'_4 + 12 * K'_5 + K'_6
    term1 = 132 * k4_prime
    term2 = 12 * k5_prime
    term3 = k6_prime
    
    total = term1 + term2 + term3

    print("The total number of successful permutations is calculated by the formula:")
    print("Total = 132 * K'_4 + 12 * K'_5 + K'_6")
    print(f"where K'_4 = {k4_prime}, K'_5 = {k5_prime}, and K'_6 = {k6_prime}.")
    print("\nCalculating each term:")
    print(f"132 * {k4_prime} = {term1}")
    print(f"12 * {k5_prime} = {term2}")
    print(f"The last term K'_6 is {term3}")
    print("\nThe sum is:")
    print(f"{term1} + {term2} + {term3} = {total}")
    print(f"\nOut of 2,985,984 possible scenarios, {total} result in the cube returning to its original configuration during the final 3 moves.")


calculate_permutations()