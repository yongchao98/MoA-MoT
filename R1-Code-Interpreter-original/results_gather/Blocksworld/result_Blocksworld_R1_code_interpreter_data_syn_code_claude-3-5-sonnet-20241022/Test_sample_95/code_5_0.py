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
    make_move(stacks, 4, 2, moves)
    
    # Step 2: Move B from stack 4 to stack 1
    make_move(stacks, 4, 1, moves)
    
    # Step 3: Move F from stack 4 to stack 2
    make_move(stacks, 4, 2, moves)
    
    # Step 4: Move C from stack 1 to stack 3
    make_move(stacks, 1, 3, moves)
    
    # Step 5: Move blocks from stack 3 to stack 2
    make_move(stacks, 3, 2, moves)  # Move D
    make_move(stacks, 3, 2, moves)  # Move E
    make_move(stacks, 3, 2, moves)  # Move A
    make_move(stacks, 3, 2, moves)  # Move G
    
    # Step 6: Move H to stack 1
    make_move(stacks, 3, 1, moves)
    
    # Step 7: Move J from stack 1 to stack 3
    make_move(stacks, 1, 3, moves)
    
    print_moves(moves)

solve_blocksworld()