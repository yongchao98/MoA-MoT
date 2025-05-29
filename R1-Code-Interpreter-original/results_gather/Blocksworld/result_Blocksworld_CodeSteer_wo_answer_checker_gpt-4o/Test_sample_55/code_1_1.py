from collections import deque

def blocksworld_bfs(initial_state, goal_state):
    # Helper function to generate new states
    def generate_moves(state):
        moves = []
        for i, stack in enumerate(state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, state))
                        new_state[i].pop()
                        new_state[j].append(block)
                        moves.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
        return moves

    # Initialize the queue with the initial state
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Generate all possible moves from the current state
        for new_state, move in generate_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))

# Define the initial and goal states
initial_state = (('J',), ('C', 'K', 'F', 'E', 'B', 'I'), ('G', 'A'), ('H', 'D'))
goal_state = (('F', 'G', 'H', 'I', 'J'), ('C', 'D', 'E', 'K'), ('A', 'B'))

# Find the sequence of moves
moves = blocksworld_bfs(initial_state, goal_state)

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")