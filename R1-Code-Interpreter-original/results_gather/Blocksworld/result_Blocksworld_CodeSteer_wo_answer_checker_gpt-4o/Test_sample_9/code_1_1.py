from collections import deque

def bfs_blocks_world(initial_state, goal_state):
    # Queue for BFS: each element is a tuple (current_state, moves)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, moves = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return moves

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
                            queue.append((new_state_tuple, moves + [f"Move {block} from {i+1} to {j+1}"]))

# Initial and goal states
initial_state = ((), ('E', 'B', 'I'), ('F', 'D'), ('C', 'H', 'G', 'A'))
goal_state = (('A', 'B', 'D', 'H'), ('I',), ('C', 'E', 'F', 'G'))

# Find the sequence of moves
moves = bfs_blocks_world(initial_state, goal_state)

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")