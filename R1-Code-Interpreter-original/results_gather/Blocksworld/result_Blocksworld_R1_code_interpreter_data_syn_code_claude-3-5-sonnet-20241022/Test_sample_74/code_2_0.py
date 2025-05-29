def print_solution():
    # Initial state
    stacks = [['B', 'G', 'F'], ['A', 'E', 'C'], ['D']]
    goal = [['D', 'F', 'G'], ['E'], ['A', 'B', 'C']]
    moves = []
    
    # Helper function to get top block
    def get_top(stack_idx):
        return stacks[stack_idx][-1] if stacks[stack_idx] else None
    
    # Helper function to check if move is valid and apply it
    def make_move(block, from_stack, to_stack):
        if not stacks[from_stack] or get_top(from_stack) != block:
            return False
        stacks[from_stack].pop()
        stacks[to_stack].append(block)
        moves.append(f"Move {block} from {from_stack+1} to {to_stack+1}")
        return True
    
    # Solution sequence
    # First, clear space for D F G in stack 1
    make_move('F', 0, 2)  # Move F from stack1 to stack3
    make_move('G', 0, 2)  # Move G from stack1 to stack3
    make_move('B', 0, 2)  # Move B from stack1 to stack3
    
    # Move D to its final position
    make_move('D', 2, 0)  # Move D from stack3 to stack1
    
    # Move F and G to their final positions
    make_move('G', 2, 0)  # Move G from stack3 to stack1
    make_move('F', 2, 0)  # Move F from stack3 to stack1
    
    # Arrange stack2 (E)
    make_move('C', 1, 2)  # Move C from stack2 to stack3
    
    # Move A to final position
    make_move('A', 1, 2)  # Move A from stack2 to stack3
    
    # Move B to final position
    make_move('B', 2, 2)  # Move B back to stack3
    
    print('\n'.join(moves))

print_solution()