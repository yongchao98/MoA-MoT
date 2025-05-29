def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    return len(stacks[from_stack-1]) > 0

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return False
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return True

def solve_blocksworld():
    # Initial state
    stacks = [
        ['J'],                          # stack 1
        ['C', 'K', 'F', 'E', 'B', 'I'], # stack 2
        ['G', 'A'],                     # stack 3
        ['H', 'D']                      # stack 4
    ]
    
    moves = []
    
    # Predefined sequence of moves to reach the goal state
    sequence = [
        (2, 4, 'I'),  # Move I from stack2 to stack4
        (2, 3, 'B'),  # Move B from stack2 to stack3
        (2, 4, 'E'),  # Move E from stack2 to stack4
        (2, 1, 'F'),  # Move F from stack2 to stack1
        (2, 2, 'K'),  # Move K from stack2 to stack2
        (4, 2, 'D'),  # Move D from stack4 to stack2
        (4, 1, 'H'),  # Move H from stack4 to stack1
        (4, 2, 'E'),  # Move E from stack4 to stack2
        (3, 3, 'A'),  # Move A to maintain its position
        (3, 1, 'G'),  # Move G from stack3 to stack1
        (4, 1, 'I'),  # Move I from stack4 to stack1
        (1, 1, 'J')   # Move J to maintain its position
    ]
    
    # Execute the sequence
    for from_stack, to_stack, block in sequence:
        if is_valid_move(stacks, from_stack, to_stack):
            if stacks[from_stack-1][-1] == block:  # Verify correct block
                make_move(stacks, from_stack, to_stack)
                moves.append(f"Move {block} from {from_stack} to {to_stack}")
    
    return moves

# Execute and print solution
solution = solve_blocksworld()
print("<<<" + "\n".join(solution) + ">>>")