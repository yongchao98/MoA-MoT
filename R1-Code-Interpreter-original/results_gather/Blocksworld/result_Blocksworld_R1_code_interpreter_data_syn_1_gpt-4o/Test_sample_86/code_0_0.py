# Initial and goal states
initial_state = {
    'stack1': ['B', 'C'],
    'stack2': ['A'],
    'stack3': ['F', 'D', 'E']
}

goal_state = {
    'stack1': ['A', 'C'],
    'stack2': ['F'],
    'stack3': ['B', 'D', 'E']
}

# Function to perform a move
def move_block(state, source, destination):
    block = state[source].pop()
    state[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the moves
moves = []

# Perform the moves
moves.append(move_block(initial_state, 'stack1', 'stack3'))  # Move C from stack1 to stack3
moves.append(move_block(initial_state, 'stack1', 'stack3'))  # Move B from stack1 to stack3
moves.append(move_block(initial_state, 'stack2', 'stack1'))  # Move A from stack2 to stack1
moves.append(move_block(initial_state, 'stack3', 'stack1'))  # Move C from stack3 to stack1
moves.append(move_block(initial_state, 'stack3', 'stack2'))  # Move B from stack3 to stack2

# Print the moves
for move in moves:
    print(move)