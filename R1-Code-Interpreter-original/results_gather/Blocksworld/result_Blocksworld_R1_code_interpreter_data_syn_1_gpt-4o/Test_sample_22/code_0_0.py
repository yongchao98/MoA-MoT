# Initial and goal states
initial_state = {
    'stack1': ['I', 'A', 'G'],
    'stack2': ['D', 'F', 'E'],
    'stack3': ['H', 'B', 'C']
}

goal_state = {
    'stack1': ['B', 'C', 'H'],
    'stack2': ['F', 'G', 'I'],
    'stack3': ['A', 'D', 'E']
}

# Function to print the moves
def print_moves(moves):
    for move in moves:
        print(f"Move {move[0]} from {move[1]} to {move[2]}")

# Function to perform the moves
def solve_blocksworld(initial, goal):
    moves = []
    
    # Move C from stack3 to stack1
    moves.append(('C', 'stack3', 'stack1'))
    initial['stack1'].append(initial['stack3'].pop())
    
    # Move B from stack3 to stack1
    moves.append(('B', 'stack3', 'stack1'))
    initial['stack1'].append(initial['stack3'].pop())
    
    # Move H from stack3 to stack1
    moves.append(('H', 'stack3', 'stack1'))
    initial['stack1'].append(initial['stack3'].pop())
    
    # Move G from stack1 to stack2
    moves.append(('G', 'stack1', 'stack2'))
    initial['stack2'].append(initial['stack1'].pop())
    
    # Move I from stack1 to stack2
    moves.append(('I', 'stack1', 'stack2'))
    initial['stack2'].append(initial['stack1'].pop())
    
    # Move A from stack1 to stack3
    moves.append(('A', 'stack1', 'stack3'))
    initial['stack3'].append(initial['stack1'].pop())
    
    # Move D from stack2 to stack3
    moves.append(('D', 'stack2', 'stack3'))
    initial['stack3'].append(initial['stack2'].pop())
    
    # Move F from stack2 to stack2 (no move needed, already in place)
    
    print_moves(moves)

solve_blocksworld(initial_state, goal_state)