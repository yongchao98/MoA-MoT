def print_moves(moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

def make_move(stacks, from_stack, to_stack, block, moves):
    # Verify the block is at the top of the source stack
    if stacks[from_stack-1] and stacks[from_stack-1][-1] == block:
        stacks[from_stack-1].pop()
        stacks[to_stack-1].append(block)
        moves.append(f"Move {block} from {from_stack} to {to_stack}")
        return True
    return False

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
    make_move(stacks, 4, 2, 'I', moves)  # Move I to stack 2
    make_move(stacks, 4, 1, 'B', moves)  # Move B to stack 1
    make_move(stacks, 4, 2, 'F', moves)  # Move F to stack 2
    
    # Step 2: Handle C and J from stack 1
    if stacks[0] and stacks[0][-1] == 'C':
        make_move(stacks, 1, 3, 'C', moves)  # Move C to stack 3
    
    if stacks[0] and stacks[0][-1] == 'J':
        make_move(stacks, 1, 3, 'J', moves)  # Move J to stack 3
    
    # Step 3: Move blocks from stack 3 to stack 2
    blocks_order = ['D', 'E', 'A', 'G']
    for block in blocks_order:
        if stacks[2] and stacks[2][-1] == block:
            make_move(stacks, 3, 2, block, moves)
    
    # Step 4: Move H to final position
    if stacks[2] and stacks[2][-1] == 'H':
        make_move(stacks, 3, 1, 'H', moves)  # Move H to stack 1
    
    print_moves(moves)

solve_blocksworld()