# Initial and goal states
initial_state = {
    'stack1': ['F'],
    'stack2': ['G', 'B', 'D'],
    'stack3': ['A', 'I', 'C', 'E', 'H']
}

goal_state = {
    'stack1': ['E', 'G', 'H'],
    'stack2': ['A', 'D', 'I'],
    'stack3': ['B', 'C', 'F']
}

# Function to print the moves
def print_moves(moves):
    for move in moves:
        print(f"Move {move[0]} from {move[1]} to {move[2]}")

# Function to perform the moves
def solve_blocksworld(initial, goal):
    moves = []
    stacks = initial.copy()

    # Move H from stack3 to stack1
    moves.append(('H', 'stack3', 'stack1'))
    stacks['stack1'].append(stacks['stack3'].pop())

    # Move E from stack3 to stack1
    moves.append(('E', 'stack3', 'stack1'))
    stacks['stack1'].append(stacks['stack3'].pop())

    # Move I from stack3 to stack2
    moves.append(('I', 'stack3', 'stack2'))
    stacks['stack2'].append(stacks['stack3'].pop())

    # Move A from stack3 to stack2
    moves.append(('A', 'stack3', 'stack2'))
    stacks['stack2'].append(stacks['stack3'].pop())

    # Move B from stack2 to stack3
    moves.append(('B', 'stack2', 'stack3'))
    stacks['stack3'].append(stacks['stack2'].pop())

    # Move G from stack2 to stack1
    moves.append(('G', 'stack2', 'stack1'))
    stacks['stack1'].append(stacks['stack2'].pop())

    # Move F from stack1 to stack3
    moves.append(('F', 'stack1', 'stack3'))
    stacks['stack3'].append(stacks['stack1'].pop())

    return moves

# Get the moves
moves = solve_blocksworld(initial_state, goal_state)

# Print the moves
print_moves(moves)