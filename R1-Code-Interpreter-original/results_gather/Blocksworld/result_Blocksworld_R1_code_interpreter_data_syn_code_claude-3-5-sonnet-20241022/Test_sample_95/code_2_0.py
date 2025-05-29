def print_moves(moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

def get_top_block(stack):
    return stack[-1] if stack else None

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
    
    # Step 1: Move I from stack 4 to stack 2
    moves.append(make_move(stacks, 4, 2))
    
    # Step 2: Move B from stack 4 to stack 1
    moves.append(make_move(stacks, 4, 1))
    
    # Step 3: Move C from stack 1 to stack 3
    moves.append(make_move(stacks, 1, 3))
    
    # Step 4: Move J from stack 1 to temporary position (stack 2)
    moves.append(make_move(stacks, 1, 2))
    
    # Step 5: Move blocks from stack 3 to stack 2 in sequence
    while stacks[2] and get_top_block(stacks[2]) in ['D', 'E', 'A', 'G']:
        moves.append(make_move(stacks, 3, 2))
    
    # Step 6: Move F from stack 4 to stack 2
    moves.append(make_move(stacks, 4, 2))
    
    # Step 7: Move H from stack 3 to stack 1
    moves.append(make_move(stacks, 3, 1))
    
    # Step 8: Move J from stack 2 to stack 3
    for i, stack in enumerate(stacks):
        if stack and 'J' in stack and stack[-1] == 'J':
            moves.append(make_move(stacks, i+1, 3))
            break
    
    print_moves(moves)

solve_blocksworld()