def print_solution():
    # Initial state
    stacks = [
        ['H', 'D'],
        ['E', 'J', 'A', 'C', 'F'],
        ['B', 'I', 'G', 'K']
    ]
    
    # Goal state
    goal = [
        ['B', 'E', 'I', 'J', 'K'],
        ['A', 'C', 'D', 'F'],
        ['G', 'H']
    ]
    
    moves = []
    
    # Helper function to find block location
    def find_block(block):
        for i, stack in enumerate(stacks):
            if block in stack and stack[-1] == block:  # Only if it's on top
                return i
        return None
    
    # Helper function to move a block
    def move_block(block, from_stack, to_stack):
        if from_stack is not None and stacks[from_stack][-1] == block:
            stacks[from_stack].pop()
            stacks[to_stack].append(block)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
            return True
        return False

    # First, clear space for goal stacks by moving blocks to stack 3 (temporary storage)
    for block in ['F', 'C', 'A', 'J']:
        loc = find_block(block)
        if loc is not None:
            move_block(block, loc, 2)

    # Move blocks to their final positions in stack 1
    for block in ['B', 'E', 'I', 'J', 'K']:
        loc = find_block(block)
        if loc is not None:
            move_block(block, loc, 0)

    # Move blocks to their final positions in stack 2
    for block in ['A', 'C', 'D', 'F']:
        loc = find_block(block)
        if loc is not None:
            move_block(block, loc, 1)

    # Move blocks to their final positions in stack 3
    for block in ['G', 'H']:
        loc = find_block(block)
        if loc is not None:
            move_block(block, loc, 2)

    print("<<<")
    print('\n'.join(moves))
    print(">>>")

print_solution()