def print_solution():
    # Initial state
    stacks = [['A', 'E'], ['D'], ['C', 'B']]
    moves = []
    
    # Helper function to get top block of a stack
    def get_top(stack_idx):
        return stacks[stack_idx][-1] if stacks[stack_idx] else None
    
    # Helper function to make a move
    def make_move(from_stack, to_stack):
        block = stacks[from_stack].pop()
        stacks[to_stack].append(block)
        moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
    
    # Step 1: Move B from stack 3 to stack 1
    if get_top(2) == 'B':  # stack3 is index 2
        make_move(2, 0)    # move to stack1 (index 0)
    
    # Step 2: Move C from stack 3 to stack 2
    if get_top(2) == 'C':  # stack3 is index 2
        make_move(2, 1)    # move to stack2 (index 1)
    
    # Step 3: Move D from stack 2 to stack 3
    if get_top(1) == 'D':  # stack2 is index 1
        make_move(1, 2)    # move to stack3 (index 2)
    
    # Print solution
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()