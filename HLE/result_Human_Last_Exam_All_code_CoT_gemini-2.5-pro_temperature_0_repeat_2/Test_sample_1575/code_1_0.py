def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    N = 100  # Total number of elements
    M = 5    # The modulus for free swaps (swap i and i+5)

    # Number of elements in each group
    num_per_group = N // M

    print(f"The {N} elements are divided into {M} groups based on their index modulo {M}.")
    print(f"Each group contains {num_per_group} elements.\n")

    # The groups are arranged in a cycle: 0-1-2-3-4-0
    # An element at index i (group k = i mod M) must move to index (N-1)-i.
    # The target group is ((N-1) - k) mod M.
    # For N=100, M=5, target group is (99-k) mod 5 = (4-k) mod 5.

    total_moves = 0
    group_moves = []
    equation_parts = []

    print("Calculating moves for each group:")
    for k in range(M):
        # Target group for elements starting in group k
        target_k = (N - 1 - k) % M

        # Calculate the shortest distance on the cycle of M groups
        dist = abs(target_k - k)
        cyclic_dist = M - dist
        shortest_dist = min(dist, cyclic_dist)

        # Moves for this group
        moves = num_per_group * shortest_dist
        group_moves.append(moves)
        equation_parts.append(str(moves))
        total_moves += moves

        print(f"  Group {k} (20 elements) -> Group {target_k}:")
        print(f"    Distance on cycle = {shortest_dist}")
        print(f"    Moves = {num_per_group} * {shortest_dist} = {moves}")

    # Print the final equation and result
    equation = " + ".join(equation_parts)
    print(f"\nTotal minimum moves = {equation} = {total_moves}")

    return total_moves

# Run the solver
final_answer = solve_reversal_moves()
print(f"\n<<<The final answer is {final_answer}>>>")
