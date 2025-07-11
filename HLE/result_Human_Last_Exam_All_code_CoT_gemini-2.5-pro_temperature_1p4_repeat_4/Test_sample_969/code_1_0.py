def solve_sequence():
    """
    Solves the sequence puzzle by generating the next four elements based on the identified pattern.
    """
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The logic to find the next elements follows a cyclical pattern.
    # The first cycle has ended with the number 2.
    # The new cycle starts with a countdown from this last number.
    
    next_elements = []

    # Start of the second cycle
    # G5: Countdown from the last element of the previous sequence (2).
    g5_start = initial_sequence[-1]
    g5 = list(range(g5_start, 0, -1))  # Countdown from 2 -> [2, 1]
    next_elements.extend(g5)
    
    # G6: Countup from the last element of G5 + 1.
    # The last element of G5 is 1. So, we count up from 1 + 1 = 2.
    # The length of the countup group is its starting number, which is 2.
    g6_start = g5[-1] + 1
    g6_length = g6_start
    g6 = list(range(g6_start, g6_start + g6_length)) # Countup from 2 for 2 steps -> [2, 3]
    
    # We only need enough elements from G6 to get a total of 4 new elements.
    needed = 4 - len(next_elements)
    next_elements.extend(g6[:needed])

    print("The original sequence is:")
    print(' '.join(map(str, initial_sequence)))
    
    print("\nFollowing the pattern, the next cycle begins:")
    print("1. A countdown from 2 gives the first two new elements: " + ' '.join(map(str, g5)))
    print("2. A countup from 2 gives the next two new elements: " + ' '.join(map(str, g6)))

    print("\nTherefore, the complete sequence with the next 4 elements is:")
    final_sequence_str = ' '.join(map(str, initial_sequence + next_elements))
    print(final_sequence_str)

solve_sequence()