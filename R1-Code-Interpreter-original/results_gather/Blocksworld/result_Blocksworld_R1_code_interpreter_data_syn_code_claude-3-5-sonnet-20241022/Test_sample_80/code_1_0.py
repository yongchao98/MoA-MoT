def print_moves(initial, goal):
    # Convert stacks to more manageable format
    current = [list(stack) for stack in initial]
    target = [list(stack) for stack in goal]
    moves = []
    
    # Helper function to find block location
    def find_block(block):
        for i, stack in enumerate(current):
            if stack and stack[-1] == block:
                return i
        return None

    # Helper function to find where a block should go
    def find_target_stack(block):
        for i, stack in enumerate(target):
            if block in stack:
                return i
        return None

    # Helper function to make a move
    def make_move(from_stack, to_stack):
        if current[from_stack]:
            block = current[from_stack].pop()
            current[to_stack].append(block)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")

    # Process each block in reverse order of goal stacks
    for stack_idx, stack in enumerate(target):
        for block in reversed(stack):
            current_loc = find_block(block)
            if current_loc is not None and current_loc != stack_idx:
                make_move(current_loc, stack_idx)

# Initial and goal states
initial = [
    ['E', 'C', 'J'],
    ['H', 'D'],
    ['G', 'I', 'A', 'B', 'F']
]

goal = [
    ['A', 'C', 'D', 'E', 'F', 'G'],
    ['B', 'H', 'J'],
    ['I']
]

print_moves(initial, goal)