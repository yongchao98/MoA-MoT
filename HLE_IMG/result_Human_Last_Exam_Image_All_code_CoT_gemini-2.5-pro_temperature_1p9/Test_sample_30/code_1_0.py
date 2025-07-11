import collections

def solve_knights_puzzles():
    """
    Solves the Knights Puzzle for 5 configurations on a 4x3 board by
    performing a Breadth-First Search (BFS) for each case.
    """

    # Define board and pre-calculate all possible knight moves for each square
    rows, cols = 4, 3
    all_moves = {}
    for r in range(rows):
        for c in range(cols):
            all_moves[(r, c)] = []
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    all_moves[(r, c)].append((nr, nc))

    # Define the 5 initial configurations from the image
    configs = {
        'A': {
            'W': frozenset([(0, 2), (1, 2), (2, 2), (3, 2)]),
            'B': frozenset([(0, 0), (1, 0), (2, 0), (3, 0)])
        },
        'B': {
            'W': frozenset([(1, 1), (3, 0), (3, 2)]),
            'B': frozenset([(0, 1), (1, 0), (2, 2)])
        },
        'C': {
            'W': frozenset([(0, 0), (2, 1)]),
            'B': frozenset([(0, 2), (1, 2)])
        },
        'D': {
            'W': frozenset([(0, 1), (2, 2)]),
            'B': frozenset([(1, 1), (3, 0)])
        },
        'E': {
            'W': frozenset([(0, 1), (0, 2), (1, 2)]),
            'B': frozenset([(0, 0), (1, 0), (1, 1)])
        }
    }

    solvable_configs = []
    print("Analyzing configurations...")

    for name, config in configs.items():
        initial_W = config['W']
        initial_B = config['B']
        goal_W = initial_B
        goal_B = initial_W

        initial_state = (initial_W, initial_B, 'W')
        queue = collections.deque([initial_state])
        visited = {initial_state}
        
        is_solvable = False
        
        # BFS explores the state space.
        # Limit search size to prevent excessive runtimes, though the state space is manageable.
        max_states_to_check = 100000 
        count = 0

        while queue and count <= max_states_to_check:
            count += 1
            current_W, current_B, turn = queue.popleft()
            
            if current_W == goal_W and current_B == goal_B:
                is_solvable = True
                break

            occupied_squares = current_W.union(current_B)
            
            active_knights = current_W if turn == 'W' else current_B
            for r_k, c_k in active_knights:
                for next_pos in all_moves[(r_k, c_k)]:
                    if next_pos not in occupied_squares:
                        if turn == 'W':
                            new_W = (current_W - {(r_k, c_k)}).union({next_pos})
                            next_state = (new_W, current_B, 'B')
                        else: # turn == 'B'
                            new_B = (current_B - {(r_k, c_k)}).union({next_pos})
                            next_state = (current_W, new_B, 'W')
                        
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append(next_state)

        if is_solvable:
            solvable_configs.append(name)
        
        print(f"Configuration {name}: {'Solvable' if is_solvable else 'Unsolvable'}")

    print("\nConclusion:")
    print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")

if __name__ == '__main__':
    solve_knights_puzzles()