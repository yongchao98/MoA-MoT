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
    
    # 1. First move C to its final position in stack 2
    moves.append(make_move(current_state, 4, 2))  # Move C from 4 to 2
    
    # 2. Move K to its final position in stack 3
    moves.append(make_move(current_state, 4, 3))  # Move K from 4 to 3
    
    # 3. Move D to its final position in stack 3
    moves.append(make_move(current_state, 4, 3))  # Move D from 4 to 3
    
    # 4. Move F to its final position in stack 2
    moves.append(make_move(current_state, 4, 2))  # Move F from 4 to 2
    
    # 5. Move blocks from stack 1
    moves.append(make_move(current_state, 1, 3))  # Move A from 1 to 3
    moves.append(make_move(current_state, 1, 3))  # Move E from 1 to 3
    moves.append(make_move(current_state, 1, 2))  # Move B from 1 to 2
    moves.append(make_move(current_state, 1, 2))  # Move I from 1 to 2
    
    # 6. Move G to its final position
    moves.append(make_move(current_state, 3, 2))  # Move G from 3 to 2
    
    # 7. Move J to its final position
    moves.append(make_move(current_state, 2, 3))  # Move J from 2 to 3
    
    # 8. Rearrange stack 2 to get correct order (B C F G I)
    moves.append(make_move(current_state, 2, 4))  # Move I to 4
    moves.append(make_move(current_state, 2, 4))  # Move G to 4
    moves.append(make_move(current_state, 2, 4))  # Move F to 4
    moves.append(make_move(current_state, 2, 3))  # Move C to 3 (temporary)
    moves.append(make_move(current_state, 2, 3))  # Move B to 3 (temporary)
    moves.append(make_move(current_state, 3, 2))  # Move B back to 2
    moves.append(make_move(current_state, 3, 2))  # Move C back to 2
    moves.append(make_move(current_state, 4, 2))  # Move F back to 2
    moves.append(make_move(current_state, 4, 2))  # Move G back to 2
    moves.append(make_move(current_state, 4, 2))  # Move I back to 2

    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()