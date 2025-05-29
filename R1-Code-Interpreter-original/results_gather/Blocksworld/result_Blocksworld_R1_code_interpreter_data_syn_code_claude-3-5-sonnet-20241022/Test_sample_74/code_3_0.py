def print_solution():
    # Initial state
    state = {
        1: ['B', 'G', 'F'],
        2: ['A', 'E', 'C'],
        3: ['D']
    }
    
    moves = []
    
    def get_top(stack_num):
        return state[stack_num][-1] if state[stack_num] else None
    
    def move_block(block, from_stack, to_stack):
        if not state[from_stack] or state[from_stack][-1] != block:
            return False
        if get_top(from_stack) == block:
            state[from_stack].pop()
            state[to_stack].append(block)
            moves.append(f"Move {block} from {from_stack} to {to_stack}")
            return True
        return False

    # Sequence of moves to reach goal state
    # First clear stack1 for D F G
    move_block('F', 1, 3)  # Move F to make space
    move_block('G', 1, 2)  # Move G temporarily to stack2
    move_block('B', 1, 3)  # Move B to its eventual destination stack3
    
    # Move C to its final destination
    move_block('C', 2, 3)  # Move C to stack3
    
    # Move A to its final destination
    move_block('A', 2, 3)  # Move A to stack3
    
    # Get D to stack1
    move_block('D', 3, 1)  # Move D to stack1
    
    # Get F and G to stack1
    move_block('F', 3, 1)  # Move F back to stack1
    move_block('G', 2, 1)  # Move G to stack1
    
    # E is already in correct position in stack2
    
    print('\n'.join(moves))

print_solution()