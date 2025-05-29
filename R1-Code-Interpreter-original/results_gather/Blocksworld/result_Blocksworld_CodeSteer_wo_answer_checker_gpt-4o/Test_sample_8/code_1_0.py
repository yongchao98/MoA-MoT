from collections import deque

def blocksworld_bfs(initial_state, goal_state):
    # Initialize the queue with the initial state and an empty move list
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, moves = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return moves

        # Generate all possible valid moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(3):  # Try to move it to each of the three stacks
                    if i != j:  # Don't move to the same stack
                        new_state = list(map(list, current_state))
                        new_state[i].pop()  # Remove the block from the current stack
                        new_state[j].append(block)  # Add the block to the new stack
                        new_state_tuple = tuple(map(tuple, new_state))

                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            queue.append((new_state_tuple, moves + [f"Move {block} from {i+1} to {j+1}"]))

# Define the initial and goal states
initial_state = (('D',), ('G', 'C', 'A', 'F'), ('B', 'I', 'E', 'H'))
goal_state = (('A', 'F', 'H'), ('C', 'G', 'I'), ('B', 'D', 'E'))

# Find the sequence of moves
moves = blocksworld_bfs(initial_state, goal_state)

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")