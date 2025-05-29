def get_moves():
    # Initial state
    stacks = [
        ['I'],
        ['G', 'C', 'H', 'A', 'E'],
        ['F', 'B'],
        ['D', 'K', 'J']
    ]
    
    # Track moves
    moves = []
    
    # Helper function to get top block
    def get_top(stack_idx):
        return stacks[stack_idx][-1] if stacks[stack_idx] else None
    
    # Helper function to make a move
    def make_move(from_idx, to_idx):
        if stacks[from_idx]:
            block = stacks[from_idx].pop()
            stacks[to_idx].append(block)
            moves.append(f"Move {block} from {from_idx + 1} to {to_idx + 1}")
    
    # First clear stack 4 as it's empty in goal state
    while stacks[3]:
        block = get_top(3)
        if block == 'J':  # J goes to stack 1
            make_move(3, 0)
        elif block == 'K':  # K goes to stack 3
            make_move(3, 2)
        elif block == 'D':  # D goes to stack 1
            make_move(3, 0)
    
    # Move E from stack 2 to stack 3
    if 'E' in stacks[1]:
        while get_top(1) != 'E':
            temp = get_top(1)
            make_move(1, 3)  # Use stack 4 as temporary storage
        make_move(1, 2)  # Move E to stack 3
    
    # Move A to stack 1
    if 'A' in stacks[1]:
        while get_top(1) != 'A':
            temp = get_top(1)
            if temp == 'H':  # H goes to stack 2
                make_move(1, 1)
            else:
                make_move(1, 3)  # Use stack 4 as temporary
        make_move(1, 0)
    
    # Move remaining blocks to their final positions
    while stacks[1] or stacks[2] or stacks[3]:
        for i in range(4):
            if not stacks[i]:
                continue
            block = get_top(i)
            if block == 'F':
                make_move(i, 0)  # F goes to stack 1
            elif block == 'B':
                make_move(i, 2)  # B goes to stack 3
            elif block == 'C':
                make_move(i, 2)  # C goes to stack 3
            elif block == 'G':
                make_move(i, 1)  # G goes to stack 2
            elif block == 'H':
                make_move(i, 1)  # H goes to stack 2
    
    return moves

# Print solution
print('\n'.join(get_moves()))