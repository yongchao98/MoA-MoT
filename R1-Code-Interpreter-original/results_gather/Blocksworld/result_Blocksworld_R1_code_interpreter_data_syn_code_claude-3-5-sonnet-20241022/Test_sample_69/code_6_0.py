def make_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # Check if source stack is empty
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")
    print()

def solve_blocksworld():
    # Initial state
    current_state = [
        ['H', 'I', 'B', 'E', 'A'],
        ['J'],
        ['G'],
        ['F', 'D', 'K', 'C']
    ]
    
    moves = []
    print("Initial state:")
    print_state(current_state)
    
    # Step by step solution
    # 1. First clear stack 4 and place blocks in their final positions
    moves.append(make_move(current_state, 4, 2))  # Move C from 4 to 2
    moves.append(make_move(current_state, 4, 3))  # Move K from 4 to 3 (final position)
    moves.append(make_move(current_state, 4, 3))  # Move D from 4 to 3 (final position)
    moves.append(make_move(current_state, 4, 2))  # Move F from 4 to 2 (final position)
    
    # 2. Move blocks from stack 1 to their final positions
    moves.append(make_move(current_state, 1, 3))  # Move A from 1 to 3 (final position)
    moves.append(make_move(current_state, 1, 3))  # Move E from 1 to 3 (final position)
    moves.append(make_move(current_state, 1, 2))  # Move B from 1 to 2 (final position)
    moves.append(make_move(current_state, 1, 2))  # Move I from 1 to 2 (final position)
    # Now only H remains in stack 1 (already in final position)
    
    # 3. Move G to its final position
    moves.append(make_move(current_state, 3, 2))  # Move G from 3 to 2 (final position)
    
    # 4. Move J to its final position
    moves.append(make_move(current_state, 2, 3))  # Move J from 2 to 3 (final position)
    
    # 5. Rearrange stack 2 if needed to get B C F G I order
    moves.append(make_move(current_state, 2, 4))  # Move I to 4 (temporary)
    moves.append(make_move(current_state, 2, 4))  # Move G to 4 (temporary)
    moves.append(make_move(current_state, 2, 4))  # Move F to 4 (temporary)
    moves.append(make_move(current_state, 2, 4))  # Move C to 4 (temporary)
    # Now B is at the bottom of stack 2 (correct position)
    moves.append(make_move(current_state, 4, 2))  # Move C back to 2
    moves.append(make_move(current_state, 4, 2))  # Move F back to 2
    moves.append(make_move(current_state, 4, 2))  # Move G back to 2
    moves.append(make_move(current_state, 4, 2))  # Move I back to 2

    print("Final state:")
    print_state(current_state)
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()