def solve_blocks():
    # Initialize states
    stacks = {
        1: [],
        2: ['B', 'A', 'C'],
        3: ['E', 'D']
    }
    
    moves = []
    
    def get_top(stack_num):
        return stacks[stack_num][-1] if stacks[stack_num] else None
    
    def move_block(from_stack, to_stack):
        if not stacks[from_stack]:
            return False
        block = stacks[from_stack].pop()
        stacks[to_stack].append(block)
        moves.append(f"Move {block} from {from_stack} to {to_stack}")
        return True

    # Solution sequence
    # Move C from stack 2 to free up B
    move_block(2, 3)
    
    # Move A from stack 2 to stack 3
    move_block(2, 3)
    
    # Move B from stack 2 to stack 1
    move_block(2, 1)
    
    # Move D from stack 3 to stack 2
    move_block(3, 2)
    
    # Move E from stack 3 to stack 1
    move_block(3, 1)
    
    # Move A from stack 3 to stack 3 (temporary)
    move_block(3, 3)
    
    # Move C from stack 3 to stack 2
    move_block(3, 2)
    
    print("\n".join(moves))

solve_blocks()