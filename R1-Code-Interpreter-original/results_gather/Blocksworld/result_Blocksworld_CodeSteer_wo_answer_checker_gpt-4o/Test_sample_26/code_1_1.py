from collections import deque

def blocks_world_bfs(initial_state, goal_state):
    # Queue for BFS: stores tuples of (current_state, moves)
    queue = deque([(initial_state, [])])
    visited = set()  # To keep track of visited states

    while queue:
        current_state, moves = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return moves

        # Mark the current state as visited
        visited.add(current_state)

        # Generate all possible moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, current_state))  # Convert to list of lists
                        new_state[i].pop()  # Remove the block from the current stack
                        new_state[j].append(block)  # Add the block to the new stack
                        new_state_tuple = tuple(map(tuple, new_state))  # Convert back to tuple of tuples

                        if new_state_tuple not in visited:
                            new_moves = moves + [f"Move {block} from stack{i+1} to stack{j+1}"]
                            queue.append((new_state_tuple, new_moves))

# Initial and goal states
initial_state = (('B', 'K', 'J', 'C'), ('F', 'H'), ('D', 'E', 'G', 'I', 'A'))
goal_state = (('C', 'D', 'F', 'I', 'K'), ('B', 'H', 'J'), ('A', 'E', 'G'))

# Find the sequence of moves
solution_moves = blocks_world_bfs(initial_state, goal_state)

# Print the solution
print("<<<" + "\n".join(solution_moves) + ">>>")