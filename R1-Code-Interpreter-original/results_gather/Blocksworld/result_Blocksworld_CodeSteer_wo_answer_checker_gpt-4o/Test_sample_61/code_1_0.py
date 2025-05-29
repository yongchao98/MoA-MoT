def dfs(initial_state, goal_state):
    stack = [(initial_state, [])]  # Stack to hold states and moves
    visited = set()  # To keep track of visited states

    while stack:
        current_state, moves = stack.pop()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return moves

        # Mark the current state as visited
        visited.add(tuple(tuple(s) for s in current_state))

        # Generate all possible moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = [list(s) for s in current_state]
                        new_state[i].pop()  # Remove block from current stack
                        new_state[j].append(block)  # Add block to new stack

                        # Convert to tuple for immutability
                        new_state_tuple = tuple(tuple(s) for s in new_state)

                        if new_state_tuple not in visited:
                            stack.append((new_state, moves + [f"Move {block} from {i+1} to {j+1}"]))

# Define initial and goal states
initial_state = [[], ['B'], ['G', 'D', 'F', 'A', 'C'], ['E']]
goal_state = [('A', 'B', 'D', 'F'), ('C', 'E'), ('G'), ()]

# Find the sequence of moves
moves = dfs(initial_state, goal_state)

# Print the moves
print("<<<" + "\n".join(moves) + ">>>")