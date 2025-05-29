def find_solution():
    # Initialize stacks
    stacks = [
        ['G', 'I', 'D'],
        ['E', 'A', 'K'],
        ['C', 'B', 'H', 'J', 'F']
    ]
    
    goal = [
        ['A', 'D', 'F'],
        ['B', 'H', 'I', 'J', 'K'],
        ['C', 'E', 'G']
    ]
    
    moves = []
    
    def get_top(stack_idx):
        return stacks[stack_idx][-1] if stacks[stack_idx] else None
    
    def move(from_idx, to_idx):
        if not stacks[from_idx]:
            return False
        block = stacks[from_idx].pop()
        stacks[to_idx].append(block)
        moves.append(f"Move {block} from {from_idx + 1} to {to_idx + 1}")
        return True
    
    # First clear stack 3 (moving F, J, H, B to appropriate positions)
    while get_top(2) == 'F':
        move(2, 0)  # Move F to stack 1
    while get_top(2) == 'J':
        move(2, 1)  # Move J to stack 2
    while get_top(2) == 'H':
        move(2, 1)  # Move H to stack 2
    while get_top(2) == 'B':
        move(2, 1)  # Move B to stack 2
    
    # Move I to stack 2
    while get_top(0) == 'D':
        move(0, 2)  # Temporarily move D
    while get_top(0) == 'I':
        move(0, 1)  # Move I to stack 2
    
    # Move G to stack 3
    while get_top(0) == 'G':
        move(0, 2)  # Move G to stack 3
    
    # Move K to final position
    while get_top(1) == 'K':
        move(1, 1)  # Keep K in stack 2
    
    # Move A to stack 1
    while get_top(1) == 'A':
        move(1, 0)  # Move A to stack 1
    
    # Move E to stack 3
    while get_top(1) == 'E':
        move(1, 2)  # Move E to stack 3
    
    # Final adjustments for D and F
    while get_top(2) == 'D':
        move(2, 0)  # Move D back to stack 1
    while get_top(0) == 'F':
        move(0, 0)  # Keep F in stack 1
    
    return moves

# Execute and print solution
print('\n'.join(find_solution()))