def solve_hat_puzzle():
    """
    Solves the hat puzzle by applying graph theory concepts.
    """
    num_people = 12

    # The team's strategy is to create a graph where the minimum vertex cover
    # is maximized. This is achieved by creating a graph where the maximum
    # independent set is minimized.
    # The smallest possible maximum independent set for a non-empty graph is 1.
    # This corresponds to a complete graph (K12), where every person is paired
    # with every other person.
    max_independent_set_size = 1

    # According to the theorem τ(G) + α(G) = |V|, the size of the
    # minimum vertex cover (τ) is the number of vertices minus the size of the
    # maximum independent set (α).
    min_vertex_cover_size = num_people - max_independent_set_size

    # The people in the minimum vertex cover know their numbers directly.
    num_know_directly = min_vertex_cover_size

    # The people not in the cover deduce their number by elimination.
    num_by_elimination = num_people - num_know_directly

    # The total number of people who are guaranteed to know their number.
    total_guaranteed = num_know_directly + num_by_elimination

    # Final equation: N = (num_people - max_independent_set_size) + (num_people - (num_people - max_independent_set_size))
    # which simplifies to N = num_people.
    # Printing the numbers in the core calculation step as requested.
    print(f"The strategy forces the leader to reveal at least a certain number of hats.")
    print(f"This number is the size of the minimum vertex cover for the graph of pairs.")
    print(f"For the optimal strategy (a complete graph), this size is calculated as:")
    print(f"Total People - Size of Maximum Independent Set = {num_people} - {max_independent_set_size} = {min_vertex_cover_size}")
    print(f"\nSo, {min_vertex_cover_size} people know their number directly.")
    print(f"The remaining {num_by_elimination} person deduces their number by elimination.")
    print(f"Total number of people guaranteed to know their number (N) = {min_vertex_cover_size} + {num_by_elimination} = {total_guaranteed}")

solve_hat_puzzle()

<<<12>>>