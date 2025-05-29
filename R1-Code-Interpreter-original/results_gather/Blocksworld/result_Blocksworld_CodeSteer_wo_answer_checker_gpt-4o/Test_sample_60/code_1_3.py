from collections import deque

# Define the initial and goal states
initial_state = {
    'stack1': ['G', 'B', 'F', 'E', 'D', 'A'],
    'stack2': [],
    'stack3': ['C']
}

goal_state = {
    'stack1': ['C', 'G'],
    'stack2': ['A', 'D', 'E'],
    'stack3': ['B', 'F']
}

# Function to check if two states are equal
def states_equal(state1, state2):
    return all(state1[key] == state2[key] for key in state1)

# Function to generate possible moves
def generate_moves(state):
    moves = []
    for source in state:
        if state[source]:  # If the source stack is not empty
            block = state[source][-1]
            for destination in state:
                if source != destination:  # Move to a different stack
                    new_state = {k: v[:] for k, v in state.items()}
                    new_state[source].pop()
                    new_state[destination].append(block)
                    moves.append((new_state, f"Move {block} from {source} to {destination}"))
    return moves

# BFS to find the sequence of moves
def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        state_tuple = tuple(tuple(stack) for stack in current_state.values())

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if states_equal(current_state, goal_state):
            return path

        for new_state, move in generate_moves(current_state):
            queue.append((new_state, path + [move]))

    return None

# Find the sequence of moves
moves = bfs(initial_state, goal_state)

# Print the sequence of moves
if moves:
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")