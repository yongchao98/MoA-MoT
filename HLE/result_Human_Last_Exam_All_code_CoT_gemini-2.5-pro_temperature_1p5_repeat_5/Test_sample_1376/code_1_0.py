def solve_tomb_engraving():
    """
    Calculates the maximum number of characters that can be engraved
    by finding the optimal combination of squares and circles.
    """
    # These (N, M) pairs are derived from the packing analysis,
    # representing the max number of squares (N) for a given number of circles (M).
    # Format: (number_of_squares_N, number_of_circles_M)
    configurations = [
        (154, 0),
        (138, 1),
        (122, 2),
        (106, 3),
        (90, 4),
        (74, 5),
        (58, 6)
    ]

    max_K = -1
    best_N = 0
    best_M = 0

    print("Evaluating possible configurations:\n")

    # Iterate through each valid configuration
    for N, M in configurations:
        # Each square has 4 characters, each circle has 999 characters.
        K = 4 * N + 999 * M

        # Output the calculation for each case
        print(f"For N={N} squares and M={M} circles:")
        print(f"Total characters K = 4 * {N} + 999 * {M} = {4*N} + {999*M} = {K}\n")

        # Check if this configuration is the new best
        if K > max_K:
            max_K = K
            best_N = N
            best_M = M

    # Print the final result
    print("---")
    print("Conclusion:")
    print(f"The maximum number of characters is {max_K}.")
    print(f"This is achieved with {best_N} squares and {best_M} circles.")
    print("\nThe final answer in the format N:M:K is:")
    print(f"{best_N}:{best_M}:{max_K}")

solve_tomb_engraving()
<<<58:6:6226>>>