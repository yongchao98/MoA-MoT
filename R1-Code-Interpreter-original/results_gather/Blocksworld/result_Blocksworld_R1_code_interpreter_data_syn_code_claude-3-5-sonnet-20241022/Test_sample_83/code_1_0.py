def is_top_block(stacks, stack_idx, block):
    return stacks[stack_idx] and stacks[stack_idx][-1] == block

def make_move(stacks, from_stack, to_stack, block):
    # Verify the move is valid
    if not is_top_block(stacks, from_stack-1, block):
        return None
    
    # Make the move
    stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def solve_blocksworld():
    # Initial state
    stacks = [
        ['G', 'A', 'E', 'C'],  # stack 1
        ['H', 'I', 'F'],       # stack 2
        ['K'],                 # stack 3
        ['J', 'B', 'D']        # stack 4
    ]
    
    moves = []
    
    # Step 1: Move D (top of stack 4)
    moves.append(make_move(stacks, 4, 3, 'D'))
    
    # Step 2: Move B (now top of stack 4)
    moves.append(make_move(stacks, 4, 3, 'B'))
    
    # Step 3: Move J (now top of stack 4)
    moves.append(make_move(stacks, 4, 2, 'J'))
    
    # Step 4: Move C (top of stack 1)
    moves.append(make_move(stacks, 1, 3, 'C'))
    
    # Step 5: Move E (now top of stack 1)
    moves.append(make_move(stacks, 1, 3, 'E'))
    
    # Step 6: Move F (top of stack 2)
    moves.append(make_move(stacks, 2, 3, 'F'))
    
    # Step 7: Move I (now top of stack 2)
    moves.append(make_move(stacks, 2, 3, 'I'))
    
    # Step 8: Move H (now top of stack 2)
    moves.append(make_move(stacks, 2, 3, 'H'))
    
    # Step 9: Move A (now top of stack 1)
    moves.append(make_move(stacks, 1, 1, 'A'))
    
    # Step 10: Move G (now top of stack 1)
    moves.append(make_move(stacks, 1, 3, 'G'))
    
    # Step 11: Move K (from stack 3)
    moves.append(make_move(stacks, 3, 3, 'K'))
    
    # Final moves to arrange stack 1
    moves.append(make_move(stacks, 3, 1, 'C'))
    
    # Filter out None moves and print
    moves = [move for move in moves if move is not None]
    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()