def print_solution():
    # Initial state
    stacks = [['H', 'G', 'A', 'C'], ['E'], ['F'], ['B', 'D']]
    moves = []
    
    # Helper function to get top block
    def get_top(stack_idx):
        return stacks[stack_idx][-1] if stacks[stack_idx] else None
    
    # Helper function to move block
    def move(from_idx, to_idx):
        block = stacks[from_idx].pop()
        stacks[to_idx].append(block)
        moves.append(f"Move {block} from {from_idx + 1} to {to_idx + 1}")
    
    # Step 1: Clear space for A-B-E-G in stack 1
    while len(stacks[0]) > 0 and stacks[0][-1] not in ['A', 'B', 'E', 'G']:
        if stacks[0][-1] == 'C':
            move(0, 1)  # Move C to stack 2
        elif stacks[0][-1] == 'H':
            move(0, 2)  # Move H to stack 3
        else:
            move(0, 3)  # Move other blocks to stack 4
    
    # Step 2: Move blocks to their final positions
    # Move B from stack 4 to stack 1
    if 'B' in stacks[3]:
        while stacks[3][-1] != 'B':
            move(3, 2)
        move(3, 0)
    
    # Move E to stack 1
    if 'E' in stacks[1]:
        move(1, 0)
    
    # Move G to stack 1 if it's not there
    for i in range(4):
        if stacks[i] and stacks[i][-1] == 'G':
            move(i, 0)
            break
    
    # Move F to stack 2
    for i in range(4):
        if stacks[i] and stacks[i][-1] == 'F':
            move(i, 1)
            break
    
    # Move D to stack 3
    for i in range(4):
        if stacks[i] and stacks[i][-1] == 'D':
            move(i, 2)
            break
    
    print('\n'.join(moves))

print_solution()