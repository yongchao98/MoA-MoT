def print_moves(moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

def get_stack_with_block(stacks, block):
    for i, stack in enumerate(stacks):
        if block in stack:
            return i + 1
    return None

def make_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def solve_blocksworld():
    # Initial state
    stacks = [
        ['J', 'C'],           # stack 1
        [],                   # stack 2
        ['H', 'G', 'A', 'E', 'D'],  # stack 3
        ['F', 'B', 'I']      # stack 4
    ]
    
    moves = []
    
    # First move I from stack 4 to stack 2
    moves.append(make_move(stacks, 4, 2))
    
    # Move B from stack 4 to stack 1
    moves.append(make_move(stacks, 4, 1))
    
    # Move F from stack 4 to stack 2
    moves.append(make_move(stacks, 4, 2))
    
    # Move C from stack 1 to stack 3
    moves.append(make_move(stacks, 1, 3))
    
    # Move J from stack 1 to stack 3 (it will be on top of C)
    moves.append(make_move(stacks, 1, 3))
    
    # Move blocks from stack 3 to stack 2 in sequence
    blocks_to_move = ['D', 'E', 'A', 'G']
    for block in blocks_to_move:
        if stacks[2] and block in stacks[2][-1]:
            moves.append(make_move(stacks, 3, 2))
    
    # Move H from stack 3 to stack 1 (on top of B)
    moves.append(make_move(stacks, 3, 1))
    
    print_moves(moves)

solve_blocksworld()