def make_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # Check if source stack is empty
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def solve_blocksworld():
    # Initial state
    current_state = [
        ['H', 'I', 'B', 'E', 'A'],  # stack 1
        ['J'],                       # stack 2
        ['G'],                       # stack 3
        ['F', 'D', 'K', 'C']        # stack 4
    ]
    
    moves = []
    
    # 1. Clear stack 4 and place blocks in temporary positions
    moves.append(make_move(current_state, 4, 1))  # Move C from 4 to 1 (temporary)
    moves.append(make_move(current_state, 4, 3))  # Move K from 4 to 3 (final)
    moves.append(make_move(current_state, 4, 3))  # Move D from 4 to 3 (final)
    moves.append(make_move(current_state, 4, 2))  # Move F from 4 to 2 (final)
    
    # 2. Move blocks from stack 1 to final positions
    moves.append(make_move(current_state, 1, 3))  # Move A from 1 to 3 (final)
    moves.append(make_move(current_state, 1, 3))  # Move E from 1 to 3 (final)
    moves.append(make_move(current_state, 1, 2))  # Move B from 1 to 2 (final)
    moves.append(make_move(current_state, 1, 2))  # Move I from 1 to 2 (final)
    moves.append(make_move(current_state, 1, 2))  # Move C from 1 to 2 (final)
    
    # 3. Move G to final position
    moves.append(make_move(current_state, 3, 2))  # Move G from 3 to 2 (final)
    
    # 4. Move J to final position
    moves.append(make_move(current_state, 2, 3))  # Move J from 2 to 3 (final)
    
    # 5. Rearrange stack 2 to get B C F G I
    moves.append(make_move(current_state, 2, 4))  # Move I to 4
    moves.append(make_move(current_state, 2, 4))  # Move G to 4
    moves.append(make_move(current_state, 2, 4))  # Move F to 4
    moves.append(make_move(current_state, 2, 4))  # Move C to 4
    # B is now at bottom of stack 2 (correct)
    moves.append(make_move(current_state, 4, 2))  # Move C to 2
    moves.append(make_move(current_state, 4, 2))  # Move F to 2
    moves.append(make_move(current_state, 4, 2))  # Move G to 2
    moves.append(make_move(current_state, 4, 2))  # Move I to 2

    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()