def print_solution():
    # Initial state
    stacks = [['A', 'E'], ['D'], ['C', 'B']]
    moves = []
    
    # Helper function to get top block of a stack
    def get_top(stack_idx):
        return stacks[stack_idx][-1] if stacks[stack_idx] else None
    
    # Helper function to make a move
    def make_move(from_stack, to_stack):
        if stacks[from_stack]:  # Check if source stack is not empty
            block = stacks[from_stack].pop()
            stacks[to_stack].append(block)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
    
    # First move: B from stack 3 to stack 1
    if 'B' in stacks[2]:  # Check if B is in stack 3
        make_move(2, 0)
    
    # Second move: C from stack 3 to stack 2
    if 'C' in stacks[2]:  # Check if C is in stack 3
        make_move(2, 1)
    
    # Third move: D from stack 2 to stack 3
    if 'D' in stacks[1]:  # Check if D is in stack 2
        make_move(1, 2)
    
    # Print solution
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()