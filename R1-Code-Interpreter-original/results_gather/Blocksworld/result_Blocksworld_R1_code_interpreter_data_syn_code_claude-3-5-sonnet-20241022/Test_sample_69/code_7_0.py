def make_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # Check if source stack is empty
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def solve_blocksworld():
    # Initial state
    current_state = [
        ['H', 'I', 'B', 'E', 'A'],
        ['J'],
        ['G'],
        ['F', 'D', 'K', 'C']
    ]
    
    moves = []
    
    # 1. First clear stack 4 and place blocks in correct positions
    moves.append(make_move(current_state, 4, 3))  # Move C to 3 (temporary)
    moves.append(make_move(current_state, 4, 3))  # Move K to 3 (final)
    moves.append(make_move(current_state, 4, 3))  # Move D to 3 (final)
    moves.append(make_move(current_state, 4, 2))  # Move F to 2 (final)
    
    # 2. Move blocks from stack 1
    moves.append(make_move(current_state, 1, 3))  # Move A to 3 (final)
    moves.append(make_move(current_state, 1, 3))  # Move E to 3 (final)
    moves.append(make_move(current_state, 1, 2))  # Move B to 2 (final)
    moves.append(make_move(current_state, 1, 2))  # Move I to 2 (final)
    # H remains in stack 1 (already in position)
    
    # 3. Move C to its final position
    moves.append(make_move(current_state, 3, 2))  # Move C to 2 (final)
    
    # 4. Move G to its final position
    moves.append(make_move(current_state, 3, 2))  # Move G to 2 (final)
    
    # 5. Finally move J to stack 3
    moves.append(make_move(current_state, 2, 3))  # Move J to 3 (final)

    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()