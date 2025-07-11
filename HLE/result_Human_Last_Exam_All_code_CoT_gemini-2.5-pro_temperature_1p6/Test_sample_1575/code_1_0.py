def calculate_reversal_moves():
    """
    Calculates the minimum moves to reverse a 100-element sequence
    with free non-adjacent swaps (i, i+5).
    """
    n = 100
    num_groups = 5

    # Create the initial sequence of types based on position index mod 5.
    # Type `k` corresponds to positions i where i % 5 == k-1.
    initial_types = [(i % num_groups) + 1 for i in range(n)]

    # Inversions are counted for pairs (i, j) where i < j and type[i] < type[j],
    # because the target order is descending (5, 4, 3, 2, 1).
    # We will store the number of inversions caused by each type.
    inversion_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    # Count inversions
    for i in range(n):
        for j in range(i + 1, n):
            # An inversion exists if the elements are out of the target order.
            # Target order is descending, so if a smaller type is before a
            # larger type, it's an inversion.
            if initial_types[i] < initial_types[j]:
                # This inversion is attributed to the first element's type.
                inversion_counts[initial_types[i]] += 1

    # The total number of moves is the sum of all inversions.
    total_moves = sum(inversion_counts.values())

    # The problem asks to show the breakdown of the calculation.
    # Type 5 elements are never smaller than others, so their count is 0.
    moves1 = inversion_counts[1]
    moves2 = inversion_counts[2]
    moves3 = inversion_counts[3]
    moves4 = inversion_counts[4]

    print("The total number of moves is the sum of inversions contributed by each type of element:")
    print(f"Moves for type 1: {moves1}")
    print(f"Moves for type 2: {moves2}")
    print(f"Moves for type 3: {moves3}")
    print(f"Moves for type 4: {moves4}")
    print("Type 5 elements do not contribute to inversions as they are the largest type.")
    print("\nFinal calculation:")
    print(f"{moves1} + {moves2} + {moves3} + {moves4} = {total_moves}")


calculate_reversal_moves()