def print_moves(moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

def make_move(stacks, from_stack, to_stack, block, moves):
    # Only move if the block is actually on top
    if stacks[from_stack-1] and stacks[from_stack-1][-1] == block:
        block = stacks[from_stack-1].pop()
        stacks[to_stack-1].append(block)
        moves.append(f"Move {block} from {from_stack} to {to_stack}")
        return True
    return False

def get_top_block(stack):
    return stack[-1] if stack else None

def find_block_location(stacks, block):
    for i, stack in enumerate(stacks):
        if stack and stack[-1] == block:
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
    
    # Step 1: Move I to stack 2
    make_move(stacks, 4, 2, 'I', moves)
    
    # Step 2: Move B to stack 1 (final position)
    make_move(stacks, 4, 1, 'B', moves)
    
    # Step 3: Move F to stack 2
    make_move(stacks, 4, 2, 'F', moves)
    
    # Step 4: Move C to stack 3
    make_move(stacks, 1, 3, 'C', moves)
    
    # Step 5: Move J to stack 3
    make_move(stacks, 1, 3, 'J', moves)
    
    # Step 6: Move blocks from stack 3 to stack 2 in sequence
    blocks_to_move = ['D', 'E', 'A', 'G']
    for block in blocks_to_move:
        if stacks[2] and block == stacks[2][-1]:
            make_move(stacks, 3, 2, block, moves)
    
    # Step 7: Move H to stack 1
    make_move(stacks, 3, 1, 'H', moves)
    
    print_moves(moves)

# Run the solution
solve_blocksworld()