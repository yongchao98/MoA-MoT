def print_moves(initial, goal):
    # Initialize stacks
    stacks = {
        1: initial[0][:],  # stack1
        2: initial[1][:],  # stack2
        3: initial[2][:]   # stack3
    }
    moves = []
    
    def get_top(stack_num):
        return stacks[stack_num][-1] if stacks[stack_num] else None
    
    def move_block(from_stack, to_stack):
        if not stacks[from_stack]:
            return False
        block = stacks[from_stack].pop()
        stacks[to_stack].append(block)
        moves.append(f"Move {block} from {from_stack} to {to_stack}")
        return True
    
    # First clear the way for A and B to go to stack1
    while get_top(2) == 'H' or get_top(2) == 'G':
        if get_top(2) == 'H':
            move_block(2, 1)
        if get_top(2) == 'G':
            move_block(2, 3)
            
    # Move B to temporary position if needed
    if get_top(2) == 'B':
        move_block(2, 3)
        
    # Move E if it's blocking
    if get_top(2) == 'E':
        move_block(2, 3)
        
    # Clear stack1 for A
    while stacks[1]:
        move_block(1, 2)
        
    # Move A to stack1
    if 'A' in stacks[3]:
        move_block(3, 1)
        
    # Move B to stack1
    for stack in [2, 3]:
        if 'B' in stacks[stack]:
            while get_top(stack) != 'B':
                move_block(stack, 3)
            move_block(stack, 1)
            break
            
    # Move F to stack1
    for stack in [1, 2]:
        if 'F' in stacks[stack]:
            while get_top(stack) != 'F':
                move_block(stack, 3)
            move_block(stack, 1)
            break

    print("\n".join(moves))

# Initial state
initial_stacks = [
    ['C', 'F'],        # stack1
    ['E', 'B', 'G', 'H'],  # stack2
    ['A', 'D']         # stack3
]

# Goal state
goal_stacks = [
    ['A', 'B', 'F'],   # stack1
    ['H'],             # stack2
    ['C', 'D', 'E', 'G']  # stack3
]

print_moves(initial_stacks, goal_stacks)