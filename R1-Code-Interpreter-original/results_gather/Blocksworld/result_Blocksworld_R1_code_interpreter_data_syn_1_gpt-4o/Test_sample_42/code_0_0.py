# Initial and goal states
initial_state = {
    'stack1': ['A'],
    'stack2': [],
    'stack3': ['B', 'C', 'D', 'E']
}

goal_state = {
    'stack1': ['D'],
    'stack2': ['B', 'E'],
    'stack3': ['A', 'C']
}

# Function to perform the moves
def perform_moves(initial, goal):
    moves = []
    
    # Move E from stack3 to stack2
    initial['stack2'].append(initial['stack3'].pop())
    moves.append("Move E from 3 to 2")
    
    # Move D from stack3 to stack1
    initial['stack1'].append(initial['stack3'].pop())
    moves.append("Move D from 3 to 1")
    
    # Move C from stack3 to stack3 (no move needed, but for clarity, it stays)
    # No move needed for C as it stays in stack3
    
    # Move B from stack3 to stack2
    initial['stack2'].append(initial['stack3'].pop())
    moves.append("Move B from 3 to 2")
    
    # Move A from stack1 to stack3
    initial['stack3'].append(initial['stack1'].pop())
    moves.append("Move A from 1 to 3")
    
    return moves

# Get the sequence of moves
moves_sequence = perform_moves(initial_state, goal_state)

# Print the moves
for move in moves_sequence:
    print(move)