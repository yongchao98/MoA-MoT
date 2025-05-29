from collections import deque

def blocksworld_bfs(initial_state, goal_state):
    # Initialize the queue with the initial state and an empty path
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
                        new_state[i].pop()
                        new_state[j].append(block)
                        new_state_tuple = tuple(map(tuple, new_state))

                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            queue.append((new_state_tuple, path + [f"Move {block} from {i+1} to {j+1}"]))

# Define the initial and goal states
initial_state = (('H', 'I', 'B', 'E', 'A'), ('J',), ('G',), ('F', 'D', 'K', 'C'))
goal_state = (('H',), ('B', 'C', 'F', 'G', 'I'), ('A', 'D', 'E', 'J', 'K'))

# Find the sequence of moves
moves = blocksworld_bfs(initial_state, goal_state)

# Print the moves
print("<<<" + "\n".join(moves) + ">>>")