# Initial and goal states
initial_state = {
    'stack1': ['I', 'J', 'F'],
    'stack2': ['B', 'G', 'E', 'C'],
    'stack3': ['H', 'D', 'A']
}

goal_state = {
    'stack1': ['B', 'C', 'H', 'I', 'J'],
    'stack2': ['D'],
    'stack3': ['A', 'E', 'F', 'G']
}

# Function to print the moves
def print_moves(moves):
    for move in moves:
        print(f"Move {move[0]} from {move[1]} to {move[2]}")

# List to store the moves
moves = []

# Move F from stack1 to stack3
moves.append(('F', 'stack1', 'stack3'))
initial_state['stack3'].append(initial_state['stack1'].pop())

# Move J from stack1 to stack3
moves.append(('J', 'stack1', 'stack3'))
initial_state['stack3'].append(initial_state['stack1'].pop())

# Move I from stack1 to stack3
moves.append(('I', 'stack1', 'stack3'))
initial_state['stack3'].append(initial_state['stack1'].pop())

# Move C from stack2 to stack1
moves.append(('C', 'stack2', 'stack1'))
initial_state['stack1'].append(initial_state['stack2'].pop())

# Move E from stack2 to stack3
moves.append(('E', 'stack2', 'stack3'))
initial_state['stack3'].append(initial_state['stack2'].pop())

# Move G from stack2 to stack3
moves.append(('G', 'stack2', 'stack3'))
initial_state['stack3'].append(initial_state['stack2'].pop())

# Move B from stack2 to stack1
moves.append(('B', 'stack2', 'stack1'))
initial_state['stack1'].append(initial_state['stack2'].pop())

# Move A from stack3 to stack2
moves.append(('A', 'stack3', 'stack2'))
initial_state['stack2'].append(initial_state['stack3'].pop())

# Move D from stack3 to stack2
moves.append(('D', 'stack3', 'stack2'))
initial_state['stack2'].append(initial_state['stack3'].pop())

# Move H from stack3 to stack1
moves.append(('H', 'stack3', 'stack1'))
initial_state['stack1'].append(initial_state['stack3'].pop())

# Move I from stack3 to stack1
moves.append(('I', 'stack3', 'stack1'))
initial_state['stack1'].append(initial_state['stack3'].pop())

# Move J from stack3 to stack1
moves.append(('J', 'stack3', 'stack1'))
initial_state['stack1'].append(initial_state['stack3'].pop())

# Print the moves
print_moves(moves)