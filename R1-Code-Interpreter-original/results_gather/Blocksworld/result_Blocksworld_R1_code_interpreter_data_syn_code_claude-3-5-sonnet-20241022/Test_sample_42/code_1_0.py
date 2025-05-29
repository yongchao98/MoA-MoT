def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty and not moving to same stack
    if not stacks[from_stack] or from_stack == to_stack:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    block = stacks[from_stack].pop()
    stacks[to_stack].append(block)
    return f"Move {block} from {from_stack + 1} to {to_stack + 1}"

def solve_blocksworld():
    # Initial state
    stacks = [['A'], [], ['B', 'C', 'D', 'E']]
    moves = []
    
    # Corrected sequence of moves
    moves.append(make_move(stacks, 2, 1))  # Move E to stack2
    moves.append(make_move(stacks, 2, 0))  # Move D to stack1
    moves.append(make_move(stacks, 2, 0))  # Move C to stack1
    moves.append(make_move(stacks, 2, 1))  # Move B to stack2
    moves.append(make_move(stacks, 0, 2))  # Move C to stack3
    moves.append(make_move(stacks, 0, 2))  # Move A to stack3
    
    print("<<<" + "\n".join(moves) + ">>>")

solve_blocksworld()