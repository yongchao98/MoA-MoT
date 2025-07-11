import collections

def solve_knights_puzzle():
    """
    Analyzes five configurations of the Knights Puzzle on a 4x3 board
    and determines which ones are solvable.
    """

    # The 4x3 board is indexed as follows:
    #  0  1  2
    #  3  4  5
    #  6  7  8
    #  9 10 11
    # Adjacency list representing all possible knight moves.
    ADJ = {
        0: [5, 7], 1: [6, 8], 2: [3, 7], 3: [2, 8, 10], 4: [9, 11],
        5: [0, 6, 10], 6: [1, 5, 11], 7: [0, 2], 8: [1, 3, 9],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # The five initial configurations from the image.
    # 'W' for white knights, 'B' for black knights.
    # Positions are represented using frozenset for efficient hashing.
    configs = {
        'A': {'W': frozenset([2, 5, 8, 11]), 'B': frozenset([0, 3, 6, 9])},
        'B': {'W': frozenset([4, 9, 11]),   'B': frozenset([1, 6, 8])},
        'C': {'W': frozenset([0, 7]),       'B': frozenset([2, 5])},
        'D': {'W': frozenset([0, 6]),       'B': frozenset([4, 10])},
        'E': {'W': frozenset([1, 2, 5]),   'B': frozenset([0, 3, 4])}
    }

    solvable_configs = []

    for name, pos_dict in configs.items():
        if is_solvable(pos_dict['W'], pos_dict['B'], ADJ):
            solvable_configs.append(name)

    print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")


def is_solvable(initial_white, initial_black, adj_list):
    """
    Determines if a puzzle configuration is solvable using Breadth-First Search.

    Args:
        initial_white (frozenset): The initial positions of white knights.
        initial_black (frozenset): The initial positions of black knights.
        adj_list (dict): The adjacency list for knight moves.

    Returns:
        bool: True if solvable, False otherwise.
    """
    # The goal state is where the positions are swapped.
    goal_white = initial_black
    goal_black = initial_white

    # The state is (current_white_positions, current_black_positions, turn).
    # White moves first.
    q = collections.deque([(initial_white, initial_black, 'W')])
    visited = {(initial_white, initial_black, 'W')}

    # BFS main loop
    while q:
        white_pos, black_pos, turn = q.popleft()

        # Check if the goal state has been reached.
        if white_pos == goal_white and black_pos == goal_black:
            return True

        if turn == 'W':
            movable_knights = white_pos
            next_turn = 'B'
        else:  # turn == 'B'
            movable_knights = black_pos
            next_turn = 'W'

        occupied_squares = white_pos.union(black_pos)

        # Generate next possible states.
        for start_pos in movable_knights:
            for end_pos in adj_list[start_pos]:
                if end_pos not in occupied_squares:
                    # Create the new set of positions for the moved knight.
                    new_movable_set = (movable_knights - {start_pos}).union({end_pos})

                    if turn == 'W':
                        next_state = (new_movable_set, black_pos, next_turn)
                    else:
                        next_state = (white_pos, new_movable_set, next_turn)

                    if next_state not in visited:
                        visited.add(next_state)
                        q.append(next_state)

    # If the queue is empty and the goal was not found, it's unsolvable.
    return False

if __name__ == "__main__":
    solve_knights_puzzle()