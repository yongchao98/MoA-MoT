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
    
    # Clear space in stack 4 (will be used as temporary storage)
    moves.append(make_move(current_state, 4, 1))  # Move C from 4 to 1 (temporary)
    moves.append(make_move(current_state, 4, 3))  # Move K from 4 to 3
    moves.append(make_move(current_state, 4, 3))  # Move D from 4 to 3
    moves.append(make_move(current_state, 4, 2))  # Move F from 4 to 2
    
    # Move blocks from stack 1 to their final positions
    moves.append(make_move(current_state, 1, 3))  # Move A from 1 to 3
    moves.append(make_move(current_state, 1, 3))  # Move E from 1 to 3
    moves.append(make_move(current_state, 1, 2))  # Move B from 1 to 2
    moves.append(make_move(current_state, 1, 4))  # Move I to 4 (temporary)
    moves.append(make_move(current_state, 1, 2))  # Move C to 2
    
    # Move G to its final position
    moves.append(make_move(current_state, 3, 2))  # Move G from 3 to 2
    
    # Move I to its final position
    moves.append(make_move(current_state, 4, 2))  # Move I from 4 to 2
    
    # Move J to its final position
    moves.append(make_move(current_state, 2, 3))  # Move J from 2 to 3

    print("Final state:")
    print_state(current_state)
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()