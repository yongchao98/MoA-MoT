def print_moves(moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

def make_move(stacks, from_stack, to_stack, moves):
    if not stacks[from_stack-1]:
        return False
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    moves.append(f"Move {block} from {from_stack} to {to_stack}")
    return True

def get_top_block(stack):
    return stack[-1] if stack else None

def find_block_stack(stacks, block):
    for i, stack in enumerate(stacks):
        if block in stack:
            return i + 1
    return None

def solve_blocksworld():
    # Initial state
    stacks = [
        ['J', 'C'],           # stack 1
        [],                   # stack 2
        ['H', 'G', 'A', 'E', 'D'],  # stack 3
        ['F', 'B', 'I']      # stack 4
    ]
    
    moves = []
    
    # Step 1: Clear stack 4 first
    # Move I to stack 2
    if get_top_block(stacks[3]) == 'I':
        make_move(stacks, 4, 2, moves)
    
    # Move B to stack 1
    if get_top_block(stacks[3]) == 'B':
        make_move(stacks, 4, 1, moves)
    
    # Move F to stack 2
    if get_top_block(stacks[3]) == 'F':
        make_move(stacks, 4, 2, moves)
    
    # Step 2: Move C and J to their final positions
    # First move C to stack 3
    if get_top_block(stacks[0]) == 'C':
        make_move(stacks, 1, 3, moves)
    
    # Then move J to stack 3
    if get_top_block(stacks[0]) == 'J':
        make_move(stacks, 1, 3, moves)
    
    # Step 3: Move blocks from stack 3 to stack 2
    blocks_order = ['D', 'E', 'A', 'G']
    for block in blocks_order:
        if get_top_block(stacks[2]) == block:
            make_move(stacks, 3, 2, moves)
    
    # Step 4: Move H to stack 1
    if get_top_block(stacks[2]) == 'H':
        make_move(stacks, 3, 1, moves)
    
    print_moves(moves)

solve_blocksworld()