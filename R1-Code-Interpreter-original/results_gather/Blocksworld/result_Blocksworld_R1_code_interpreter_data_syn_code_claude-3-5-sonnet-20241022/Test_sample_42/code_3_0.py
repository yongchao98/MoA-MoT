def make_move(stacks, from_stack, to_stack):
    if not stacks[from_stack]:
        return None
    block = stacks[from_stack].pop()
    stacks[to_stack].append(block)
    return f"Move {block} from {from_stack + 1} to {to_stack + 1}"

def solve_blocksworld():
    # Initial state
    stacks = [['A'], [], ['B', 'C', 'D', 'E']]
    moves = []
    
    # Step 1: Move E (top block) from stack3 to stack2
    moves.append(make_move(stacks, 2, 1))  # Move E from 3 to 2
    
    # Step 2: Move D from stack3 to stack1 (temporary)
    moves.append(make_move(stacks, 2, 0))  # Move D from 3 to 1
    
    # Step 3: Move C from stack3 to stack2 (temporary)
    moves.append(make_move(stacks, 2, 1))  # Move C from 3 to 2
    
    # Step 4: Move B from stack3 to stack2
    moves.append(make_move(stacks, 2, 1))  # Move B from 3 to 2
    
    # Step 5: Move C from stack2 to stack3
    moves.append(make_move(stacks, 1, 2))  # Move C from 2 to 3
    
    # Step 6: Move A from stack1 to stack3
    moves.append(make_move(stacks, 0, 2))  # Move A from 1 to 3
    
    # Step 7: Move D from stack1 to its final position in stack1 (no move needed)
    
    print("<<<" + "\n".join([move for move in moves if move is not None]) + ">>>")

solve_blocksworld()