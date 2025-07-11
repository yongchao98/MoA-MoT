import collections

def solve_sokoban():
    """
    Solves a Sokoban puzzle to find the shortest path for the boulder to the goal.

    The function uses a Breadth-First Search (BFS) to explore the state space,
    where a state is defined by the player and boulder positions. It finds all
    paths of the shortest possible length and then applies tie-breaking rules:
    1. Fewest changes in direction.
    2. Alphabetically first path.
    """
    grid_str = """
........
..T.....
........
.X......
........
.....O..
........
........
"""
    grid = [list(row) for row in grid_str.strip().split('\n')]
    GRID_SIZE = 8

    # Find initial positions of Player (T), Boulder (O), and Goal (X)
    player_pos, boulder_pos, goal_pos = None, None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == 'T':
                player_pos = (r, c)
            elif char == 'O':
                boulder_pos = (r, c)
            elif char == 'X':
                goal_pos = (r, c)

    # Initial state is (player_row, player_col, boulder_row, boulder_col)
    initial_state = (player_pos[0], player_pos[1], boulder_pos[0], boulder_pos[1])

    # BFS setup
    # The queue stores tuples of (state, path_string)
    queue = collections.deque([(initial_state, "")])
    visited = {initial_state}
    solutions = []
    min_len = float('inf')

    # Moves are ordered alphabetically ('d', 'l', 'r', 'u') for the lexicographical tie-breaker
    moves = collections.OrderedDict([
        ('d', (1, 0)),
        ('l', (0, -1)),
        ('r', (0, 1)),
        ('u', (-1, 0)),
    ])

    while queue:
        current_state, path = queue.popleft()
        (pr, pc, br, bc) = current_state

        # Pruning optimization: if a solution is found, no need to explore longer paths
        if len(path) > min_len:
            continue

        # Check if the boulder is at the goal
        if (br, bc) == goal_pos:
            if len(path) < min_len:
                min_len = len(path)
                solutions = [path]
            elif len(path) == min_len:
                solutions.append(path)
            # Continue searching to find all solutions of this shortest length
            continue

        # Explore all possible moves from the current state
        for move_char, (dr, dc) in moves.items():
            npr, npc = pr + dr, pc + dc  # new player position

            # Check if the player's move is within the grid boundaries
            if not (0 <= npr < GRID_SIZE and 0 <= npc < GRID_SIZE):
                continue

            # Case 1: Player moves into the boulder's space (a push)
            if (npr, npc) == (br, bc):
                nbr, nbc = br + dr, bc + dc  # new boulder position

                # Check if the new boulder position is valid (within the grid)
                if not (0 <= nbr < GRID_SIZE and 0 <= nbc < GRID_SIZE):
                    continue

                new_state = (npr, npc, nbr, nbc)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

            # Case 2: Player moves into an empty space
            else:
                new_state = (npr, npc, br, bc)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

    # After BFS, process the collected solutions based on tie-breaking rules
    if not solutions:
        final_solution = "No solution found."
    elif len(solutions) == 1:
        final_solution = solutions[0]
    else:
        # Tie-breaking logic for multiple shortest paths
        def count_direction_changes(p):
            if not p:
                return 0
            changes = 0
            for i in range(1, len(p)):
                if p[i] != p[i-1]:
                    changes += 1
            return changes

        # Sort solutions first by number of direction changes, then alphabetically
        solutions.sort(key=lambda p: (count_direction_changes(p), p))
        final_solution = solutions[0]
        
    print(final_solution)

solve_sokoban()