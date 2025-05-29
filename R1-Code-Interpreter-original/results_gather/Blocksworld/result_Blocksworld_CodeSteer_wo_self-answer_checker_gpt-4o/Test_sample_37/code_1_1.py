from collections import deque

def bfs_blocks_world(initial_state, goal_state):
    # Queue for BFS: stores (current_state, path_to_state)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Generate all possible moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, current_state))
                        new_state[i].pop()  # Remove block from current stack
                        new_state[j].append(block)  # Add block to new stack
                        new_state_tuple = tuple(map(tuple, new_state))

                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            queue.append((new_state_tuple, path + [f"Move {block} from {i+1} to {j+1}"]))

    return None  # If no solution is found

# Initial and goal states
initial_state = (('A', 'E', 'F'), ('B', 'C'), ('D', 'H', 'G'), ('I',))
goal_state = (('B', 'C', 'D', 'E'), ('A', 'F', 'G', 'H'), ('I',))

# Find the solution
solution = bfs_blocks_world(initial_state, goal_state)

# Print the solution
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")