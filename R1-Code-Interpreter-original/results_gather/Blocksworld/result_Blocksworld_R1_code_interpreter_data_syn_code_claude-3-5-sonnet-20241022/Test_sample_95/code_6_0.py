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
    
    # Clear stack 4 first (top to bottom: I, B, F)
    make_move(stacks, 4, 2, moves)  # Move I to stack 2
    make_move(stacks, 4, 1, moves)  # Move B to stack 1 (where it needs to be)
    make_move(stacks, 4, 2, moves)  # Move F to stack 2
    
    # Move C and J to stack 3
    make_move(stacks, 1, 3, moves)  # Move C to stack 3
    make_move(stacks, 1, 3, moves)  # Move J to stack 3 (will be on top of C)
    
    # Move blocks from stack 3 to stack 2 in sequence
    # Current stack 3 (bottom to top): H, G, A, E, D
    make_move(stacks, 3, 2, moves)  # Move D to stack 2
    make_move(stacks, 3, 2, moves)  # Move E to stack 2
    make_move(stacks, 3, 2, moves)  # Move A to stack 2
    make_move(stacks, 3, 2, moves)  # Move G to stack 2
    
    # Final move: H to stack 1
    make_move(stacks, 3, 1, moves)  # Move H to stack 1 (on top of B)
    
    print_moves(moves)

solve_blocksworld()