def print_stacks(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {stack}")

def move_block(stacks, from_stack, to_stack):
    # Convert to 0-based indexing
    from_idx = from_stack - 1
    to_idx = to_stack - 1
    
    # Check if source stack has blocks
    if not stacks[from_idx]:
        return None
    
    # Get the top block and move it
    block = stacks[from_idx].pop()
    stacks[to_idx].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def execute_solution():
    # Initialize stacks
    stacks = [
        ['A', 'G'],        # Stack 1
        ['F', 'E'],        # Stack 2
        ['C'],             # Stack 3
        ['H', 'J', 'I', 'B', 'D']  # Stack 4
    ]
    
    moves = []
    
    # Predefined sequence of moves
    sequence = [
        (4, 3),  # D from 4 to 3
        (4, 1),  # B from 4 to 1
        (4, 3),  # I from 4 to 3
        (4, 2),  # J from 4 to 2
        (4, 3),  # H from 4 to 3
        (1, 2),  # G from 1 to 2
        (2, 3),  # F from 2 to 3
        (3, 2)   # C from 3 to 2
    ]
    
    # Execute each move and store the result
    for from_stack, to_stack in sequence:
        move = move_block(stacks, from_stack, to_stack)
        if move:
            moves.append(move)
            # Debug: Print state after each move
            # print(f"\nAfter {move}:")
            # print_stacks(stacks)
    
    # Print all moves
    print('\n'.join(moves))

# Run the solution
execute_solution()