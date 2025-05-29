def print_moves(moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

def get_top_block(stack):
    return stack[-1] if stack else None

def make_move(stacks, from_stack, to_stack, moves):
    if not stacks[from_stack-1]:
        return False
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    moves.append(f"Move {block} from {from_stack} to {to_stack}")
    return True

def solve_blocksworld():
    # Initial state
    stacks = [
        ['J', 'C'],           # stack 1
        [],                   # stack 2
        ['H', 'G', 'A', 'E', 'D'],  # stack 3
        ['F', 'B', 'I']      # stack 4
    ]
    
    moves = []
    
    # Step 1: Clear stack 4
    make_move(stacks, 4, 2, moves)  # Move I to stack 2
    make_move(stacks, 4, 1, moves)  # Move B to stack 1
    make_move(stacks, 4, 2, moves)  # Move F to stack 2
    
    # Step 2: Move C to stack 3
    if get_top_block(stacks[0]) == 'C':
        make_move(stacks, 1, 3, moves)
    else:
        # If J is on top, move it temporarily to stack 2
        make_move(stacks, 1, 2, moves)  # Move J to stack 2
        make_move(stacks, 1, 3, moves)  # Move C to stack 3
        make_move(stacks, 2, 3, moves)  # Move J back to stack 3
    
    # Step 3: Move blocks from stack 3 to stack 2
    blocks_to_move = ['D', 'E', 'A', 'G']
    for _ in range(len(blocks_to_move)):
        if stacks[2] and get_top_block(stacks[2]) in blocks_to_move:
            make_move(stacks, 3, 2, moves)
    
    # Step 4: Move H to stack 1
    if get_top_block(stacks[2]) == 'H':
        make_move(stacks, 3, 1, moves)
    
    print_moves(moves)

solve_blocksworld()