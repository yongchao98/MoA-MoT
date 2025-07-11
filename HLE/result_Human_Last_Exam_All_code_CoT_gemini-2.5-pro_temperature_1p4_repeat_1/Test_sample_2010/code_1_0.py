from collections import deque

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with tie-breaking rules.
    """
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # Initial state for the BFS queue: (player_pos, boulder_pos, path_string)
    initial_state = (player_start, boulder_start, "")
    queue = deque([initial_state])

    # Visited set to store ((player_pos), (boulder_pos)) tuples
    visited = {(player_start, boulder_start)}

    # A list to store all solutions of the shortest length found
    solutions = []
    min_path_len = -1

    # Map moves to coordinate changes
    moves = {'u': (-1, 0), 'd': (1, 0), 'l': (0, -1), 'r': (0, 1)}
    move_order = ['u', 'd', 'l', 'r'] # Alphabetical order for exploration

    while queue:
        player_pos, boulder_pos, path = queue.popleft()

        # If we have found solutions, and the current path is longer, we can stop.
        if min_path_len != -1 and len(path) >= min_path_len:
            continue

        # Check if the boulder is at the goal
        if boulder_pos == goal_pos:
            solutions.append(path)
            if min_path_len == -1:
                min_path_len = len(path)
            continue

        # Explore possible moves from the current state
        for move_char in move_order:
            dr, dc = moves[move_char]
            
            # Calculate potential new player position
            new_player_r, new_player_c = player_pos[0] + dr, player_pos[1] + dc

            # Check if the new player position is within the grid
            if not (0 <= new_player_r < grid_size and 0 <= new_player_c < grid_size):
                continue
            
            new_player_pos = (new_player_r, new_player_c)
            new_boulder_pos = boulder_pos
            
            # Check if the move is a push
            if new_player_pos == boulder_pos:
                # Calculate potential new boulder position
                new_boulder_r, new_boulder_c = boulder_pos[0] + dr, boulder_pos[1] + dc

                # Check if the new boulder position is within the grid
                if not (0 <= new_boulder_r < grid_size and 0 <= new_boulder_c < grid_size):
                    continue # Cannot push boulder out of bounds
                
                new_boulder_pos = (new_boulder_r, new_boulder_c)

            # Create the new state
            new_state = (new_player_pos, new_boulder_pos)

            if new_state not in visited:
                visited.add(new_state)
                new_path = path + move_char
                queue.append((new_player_pos, new_boulder_pos, new_path))

    # After the BFS, apply tie-breaking rules to the found solutions
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

    # Find the best path among the solutions based on tie-breakers
    best_path = ""
    min_changes = float('inf')
    
    # Sort solutions alphabetically to handle the final tie-breaker implicitly
    solutions.sort()

    for path in solutions:
        changes = count_direction_changes(path)
        if changes < min_changes:
            min_changes = changes
            best_path = path

    print(best_path)

if __name__ == "__main__":
    solve_sokoban()