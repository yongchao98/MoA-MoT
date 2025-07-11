def solve_unique_factorization_lengths():
    """
    Calculates the size of the set of specified rings that have unique factorization lengths.
    This corresponds to counting Half-Factorial Domains (HFDs) in the given set.
    """

    # Part 1: Count the rings of integers O_{Q(sqrt(-d))} that are HFDs.
    # This is equivalent to counting the imaginary quadratic fields Q(sqrt(-d))
    # with class number 1 or 2. d must be square-free.

    # d for which class number h(Q(sqrt(-d))) = 1 (9 values)
    # These are the Heegner numbers. All are square-free.
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]

    # d for which class number h(Q(sqrt(-d))) = 2 (18 values)
    # All are square-free.
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    # The set of all d for which O_{Q(sqrt(-d))} is an HFD is the union of the two lists.
    d_integrally_closed_hfd = sorted(d_h1 + d_h2)
    count_integrally_closed = len(d_integrally_closed_hfd)

    print(f"The d values for which the ring of integers O_Q(sqrt(-d)) has unique factorization lengths (class number 1 or 2) are:")
    print(d_integrally_closed_hfd)
    print(f"Number of such rings: {count_integrally_closed}\n")


    # Part 2: Count the rings Z[sqrt(-d)] that are not integrally closed and are HFDs.
    # This occurs for d > 0, square-free, and d = 3 (mod 4).
    # For such d, Z[sqrt(-d)] is an HFD if and only if:
    # a) d = 3, OR
    # b) d > 3, d = 7 (mod 8) and h(Q(sqrt(-d))) is 1 or 2.

    d_non_integrally_closed_hfd = []

    # Case a) d = 3
    # 3 is square-free and 3 = 3 (mod 4). This gives one ring.
    d_non_integrally_closed_hfd.append(3)

    # Case b) Check the d values from Part 1
    # We need d > 3 and d = 7 (mod 8)
    for d in d_integrally_closed_hfd:
        if d > 3 and d % 8 == 7:
            d_non_integrally_closed_hfd.append(d)

    d_non_integrally_closed_hfd = sorted(list(set(d_non_integrally_closed_hfd)))
    count_non_integrally_closed = len(d_non_integrally_closed_hfd)

    print(f"The d values for which the non-integrally closed ring Z[sqrt(-d)] has unique factorization lengths are:")
    print(d_non_integrally_closed_hfd)
    print(f"Number of such rings: {count_non_integrally_closed}\n")

    # Part 3: Calculate the total size.
    # The set of integrally closed rings and the set of non-integrally closed rings are disjoint.
    # So, the total size of the union is the sum of the sizes of the two sets.
    total_size = count_integrally_closed + count_non_integrally_closed
    
    # Final output as an equation
    print("The total size is the sum of the number of rings from both sets.")
    print(f"{count_integrally_closed} + {count_non_integrally_closed} = {total_size}")
    
    return total_size

# Execute the function and print the final answer in the required format.
final_answer = solve_unique_factorization_lengths()
print(f"\n<<< {final_answer} >>>")