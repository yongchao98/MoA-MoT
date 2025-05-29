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

# Function to move a block from one stack to another
def move_block(state, source, destination):
    block = state[source].pop()
    state[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the sequence of moves
moves = []

# Perform the moves to reach the goal state
moves.append(move_block(initial_state, 'stack1', 'stack2'))  # Move A
moves.append(move_block(initial_state, 'stack1', 'stack2'))  # Move D
moves.append(move_block(initial_state, 'stack1', 'stack2'))  # Move E
moves.append(move_block(initial_state, 'stack1', 'stack3'))  # Move F
moves.append(move_block(initial_state, 'stack1', 'stack3'))  # Move B
moves.append(move_block(initial_state, 'stack3', 'stack1'))  # Move C
moves.append(move_block(initial_state, 'stack1', 'stack2'))  # Move G
moves.append(move_block(initial_state, 'stack3', 'stack2'))  # Move B
moves.append(move_block(initial_state, 'stack3', 'stack2'))  # Move F

# Print the sequence of moves
print("<<<" + "\n".join(moves) + ">>>")