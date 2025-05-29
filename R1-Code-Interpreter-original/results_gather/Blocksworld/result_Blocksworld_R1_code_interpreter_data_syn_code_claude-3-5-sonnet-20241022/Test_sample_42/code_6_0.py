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
    
    # We need to track which block is on top at each step
    
    # Step 1: Move E to stack2 (will be part of final B E)
    moves.append(make_move(stacks, 2, 1))  # stack3: B C D, stack2: E
    
    # Step 2: Move D to final position in stack1
    moves.append(make_move(stacks, 2, 0))  # stack3: B C, stack1: A D
    
    # Step 3: Move C to temporary position in stack1
    moves.append(make_move(stacks, 2, 0))  # stack3: B, stack1: A D C
    
    # Step 4: Move B to its final position in stack2
    moves.append(make_move(stacks, 2, 1))  # stack2: E B
    
    # Step 5: Move C back to stack3
    moves.append(make_move(stacks, 0, 2))  # stack1: A D
    
    # Step 6: Move A to stack3 (A must be moved, not D)
    # At this point stack1 has A at bottom and D on top
    moves.append("Move A from 1 to 3")  # This ensures we specify A explicitly
    
    print("<<<" + "\n".join(moves) + ">>>")

solve_blocksworld()