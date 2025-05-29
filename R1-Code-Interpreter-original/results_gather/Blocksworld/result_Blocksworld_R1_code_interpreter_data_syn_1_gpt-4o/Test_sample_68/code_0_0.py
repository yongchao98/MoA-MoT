# Initial and goal states
initial_state = {
    'stack1': ['H', 'D'],
    'stack2': ['E', 'J', 'A', 'C', 'F'],
    'stack3': ['B', 'I', 'G', 'K']
}

goal_state = {
    'stack1': ['B', 'E', 'I', 'J', 'K'],
    'stack2': ['A', 'C', 'D', 'F'],
    'stack3': ['G', 'H']
}

# Function to print the move
def move(block, source, destination):
    print(f"Move {block} from {source} to {destination}")

# Function to perform the moves
def solve_blocks(initial, goal):
    # Move K from stack3 to stack1
    move('K', 'stack3', 'stack1')
    initial['stack3'].pop()
    initial['stack1'].append('K')
    
    # Move G from stack3 to stack3 (no move needed, already in place)
    
    # Move I from stack3 to stack1
    move('I', 'stack3', 'stack1')
    initial['stack3'].pop()
    initial['stack1'].append('I')
    
    # Move B from stack3 to stack1
    move('B', 'stack3', 'stack1')
    initial['stack3'].pop()
    initial['stack1'].append('B')
    
    # Move F from stack2 to stack3
    move('F', 'stack2', 'stack3')
    initial['stack2'].pop()
    initial['stack3'].append('F')
    
    # Move C from stack2 to stack2 (no move needed, already in place)
    
    # Move A from stack2 to stack2 (no move needed, already in place)
    
    # Move J from stack2 to stack1
    move('J', 'stack2', 'stack1')
    initial['stack2'].pop()
    initial['stack1'].append('J')
    
    # Move E from stack2 to stack1
    move('E', 'stack2', 'stack1')
    initial['stack2'].pop()
    initial['stack1'].append('E')
    
    # Move D from stack1 to stack2
    move('D', 'stack1', 'stack2')
    initial['stack1'].pop()
    initial['stack2'].append('D')
    
    # Move H from stack1 to stack3
    move('H', 'stack1', 'stack3')
    initial['stack1'].pop()
    initial['stack3'].append('H')

# Execute the solution
solve_blocks(initial_state, goal_state)