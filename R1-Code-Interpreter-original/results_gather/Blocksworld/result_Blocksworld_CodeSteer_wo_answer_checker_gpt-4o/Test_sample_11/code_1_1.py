# Define the initial and goal states
initial_state = {
    'stack1': [],
    'stack2': ['E'],
    'stack3': ['B', 'A', 'F', 'J'],
    'stack4': ['C', 'G', 'D', 'H', 'I']
}

goal_state = {
    'stack1': ['C', 'D'],
    'stack2': ['E', 'F', 'I'],
    'stack3': ['A', 'B', 'G', 'H', 'J']
}

# Function to perform a move
def move_block(state, source, destination):
    block = state[source].pop()
    state[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# Function to check if the current state matches the goal state
def is_goal_state(state, goal_state):
    return state == goal_state

# List to store the sequence of moves
moves = []

# Simulate the moves
moves.append(move_block(initial_state, 'stack4', 'stack1'))  # Move I from stack4 to stack2
moves.append(move_block(initial_state, 'stack4', 'stack2'))  # Move H from stack4 to stack3
moves.append(move_block(initial_state, 'stack4', 'stack1'))  # Move D from stack4 to stack1
moves.append(move_block(initial_state, 'stack4', 'stack3'))  # Move G from stack4 to stack3
moves.append(move_block(initial_state, 'stack4', 'stack1'))  # Move C from stack4 to stack1
moves.append(move_block(initial_state, 'stack3', 'stack2'))  # Move F from stack3 to stack2
moves.append(move_block(initial_state, 'stack2', 'stack2'))  # Move E from stack2 to stack2
moves.append(move_block(initial_state, 'stack3', 'stack3'))  # Move J from stack3 to stack3
moves.append(move_block(initial_state, 'stack3', 'stack3'))  # Move A from stack3 to stack3
moves.append(move_block(initial_state, 'stack3', 'stack3'))  # Move B from stack3 to stack3

# Check if the goal state is reached
if is_goal_state(initial_state, goal_state):
    print("Goal state reached!")
else:
    print("Goal state not reached.")

# Print the sequence of moves
print("Moves:")
for move in moves:
    print(move)