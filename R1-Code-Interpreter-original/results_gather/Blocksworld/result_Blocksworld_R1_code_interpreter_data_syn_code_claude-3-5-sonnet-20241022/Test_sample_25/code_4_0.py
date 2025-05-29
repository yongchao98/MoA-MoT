def print_state(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {' '.join(stack)}")

def is_valid_move(stacks, from_stack, to_stack):
    if from_stack < 1 or from_stack > len(stacks) or to_stack < 1 or to_stack > len(stacks):
        return False
    if not stacks[from_stack-1]:  # Empty source stack
        return False
    return True

def make_move(stacks, moves, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return False
    
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    moves.append(f"Move {block} from {from_stack} to {to_stack}")
    return True

def solve_blocksworld():
    # Initialize stacks
    stacks = [
        ['A', 'G'],        # Stack 1
        ['F', 'E'],        # Stack 2
        ['C'],             # Stack 3
        ['H', 'J', 'I', 'B', 'D']  # Stack 4
    ]
    
    moves = []
    
    # Step 1: Move D from stack 4 to 3
    make_move(stacks, moves, 4, 3)
    
    # Step 2: Move B from stack 4 to 1
    make_move(stacks, moves, 4, 1)
    
    # Step 3: Move I from stack 4 to 3
    make_move(stacks, moves, 4, 3)
    
    # Step 4: Move J from stack 4 to 2
    make_move(stacks, moves, 4, 2)
    
    # Step 5: Move H from stack 4 to 3
    make_move(stacks, moves, 4, 3)
    
    # Step 6: Move G from stack 1 to 2
    make_move(stacks, moves, 1, 2)
    
    # Step 7: Move E from stack 2 to 2 (already there)
    
    # Step 8: Move F from stack 2 to 3
    make_move(stacks, moves, 2, 3)
    
    # Step 9: Move C from stack 3 to 2
    make_move(stacks, moves, 3, 2)
    
    print('\n'.join(moves))

# Run the solution
solve_blocksworld()