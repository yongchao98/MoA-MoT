def demonstrate_maximal_rank_calculation(d):
    """
    This function demonstrates the calculation for the maximal rank of the
    complementary channel of a d-dimensional Pauli channel.

    The rank is determined by the size of the difference set S-S, where S
    is the set of indices (a,b) for which the Pauli operator U_ab has a
    non-zero probability p_ab. To maximize the rank, we choose a channel
    where all p_ab > 0, so S is the entire group G = Z_d x Z_d.

    Args:
        d (int): The dimension of the qudit system.
    """
    if not isinstance(d, int) or d < 2:
        print("Please provide an integer dimension d >= 2.")
        return

    # 1. The group G = Z_d x Z_d is represented by pairs of integers (i, j).
    # To maximize the rank, we choose S to be the entire group G.
    S = []
    for i in range(d):
        for j in range(d):
            S.append((i, j))

    # 2. Compute the difference set S-S = {s1 - s2 for s1, s2 in S}.
    # Subtraction is component-wise modulo d.
    S_minus_S = set()
    for s1 in S:
        for s2 in S:
            diff = ((s1[0] - s2[0]) % d, (s1[1] - s2[1]) % d)
            S_minus_S.add(diff)

    # 3. The maximal rank is the size of this set.
    maximal_rank = len(S_minus_S)

    # 4. Output the results, showing the final equation.
    print(f"For a qudit of dimension d = {d}:")
    print(f"The maximal rank of the complementary channel is given by d^2.")
    print("This is demonstrated by taking the set of active Pauli operators S to be the full group Z_d x Z_d.")
    print(f"The size of the resulting difference set |S-S| is {maximal_rank}.")
    print("\nThe final equation for the maximal rank is:")
    
    # Print each number in the final equation
    power = 2
    print(f"{d} ** {power} = {maximal_rank}")

# Example: Run the demonstration for a 4-dimensional system (ququart)
demonstrate_maximal_rank_calculation(4)
