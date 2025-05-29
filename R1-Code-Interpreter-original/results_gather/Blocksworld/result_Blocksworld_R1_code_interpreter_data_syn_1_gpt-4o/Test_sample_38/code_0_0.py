# Initial and goal states
initial_state = {
    'stack1': ['G', 'A', 'F', 'E', 'C'],
    'stack2': ['J', 'H', 'B'],
    'stack3': ['I', 'D']
}

goal_state = {
    'stack1': ['A', 'B', 'D', 'F', 'I'],
    'stack2': ['C', 'H'],
    'stack3': ['E', 'G', 'J']
}

# Function to perform the moves
def perform_moves(initial, goal):
    moves = []
    stacks = {1: initial['stack1'], 2: initial['stack2'], 3: initial['stack3']}
    
    # Helper function to find the stack containing a block
    def find_stack(block):
        for key, stack in stacks.items():
            if stack and stack[-1] == block:
                return key
        return None
    
    # Move blocks to achieve the goal state
    for stack_num, goal_stack in goal.items():
        for block in reversed(goal_stack):
            current_stack = find_stack(block)
            if current_stack is not None and current_stack != stack_num:
                stacks[current_stack].pop()
                stacks[stack_num].append(block)
                moves.append(f"Move {block} from {current_stack} to {stack_num}")
    
    return moves

# Get the sequence of moves
moves_sequence = perform_moves(initial_state, goal_state)

# Print the moves
print("\n".join(moves_sequence))