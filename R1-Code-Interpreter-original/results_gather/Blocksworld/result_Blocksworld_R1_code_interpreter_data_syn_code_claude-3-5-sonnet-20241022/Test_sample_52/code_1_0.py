def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # source stack is empty
        return False
    return True

def make_move(stacks, move):
    from_stack, to_stack = move
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = [stack[:] for stack in stacks]
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def solve_blocks():
    # Initial and goal states
    initial = [['I', 'J', 'F'], ['B', 'G', 'E', 'C'], ['H', 'D', 'A']]
    goal = [['B', 'C', 'H', 'I', 'J'], ['D'], ['A', 'E', 'F', 'G']]
    
    # Manual solution based on analysis
    moves = [
        (2, 3),  # Move C from 2 to 3
        (2, 3),  # Move E from 2 to 3
        (2, 3),  # Move G from 2 to 3
        (2, 1),  # Move B from 2 to 1
        (1, 3),  # Move F from 1 to 3
        (1, 1),  # Move J from 1 to 1
        (3, 2),  # Move D from 3 to 2
        (3, 1),  # Move H from 3 to 1
        (1, 1),  # Move I from 1 to 1
    ]
    
    current = initial
    solution = []
    
    for move in moves:
        from_stack, to_stack = move
        if current[from_stack-1]:  # Check if source stack is not empty
            block = current[from_stack-1][-1]  # Get top block
            new_state = make_move(current, move)
            if new_state:
                solution.append(f"Move {block} from {from_stack} to {to_stack}")
                current = new_state

    print("<<<" + "\n".join(solution) + ">>>")

solve_blocks()