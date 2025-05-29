def get_top_block(stack):
    return stack[-1] if stack else None

def make_move(stacks, from_stack, to_stack, moves):
    if not stacks[from_stack-1]:  # Check if source stack is empty
        return False
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    moves.append(f"Move {block} from {from_stack} to {to_stack}")
    return True

def solve_blocksworld(initial, goal):
    current = [list(stack) for stack in initial]
    moves = []
    
    # Step 1: Build stack 1 (E, G, H)
    # First move H temporarily
    make_move(current, 3, 2, moves)  # Move H to stack 2
    
    # Move E to final position
    make_move(current, 3, 1, moves)  # Move E to stack 1
    
    # Move blocks to get G
    make_move(current, 2, 3, moves)  # Move D
    make_move(current, 2, 3, moves)  # Move B
    make_move(current, 2, 1, moves)  # Move G to stack 1
    
    # Move H to complete stack 1
    make_move(current, 2, 1, moves)  # Move H to stack 1
    
    # Step 2: Build stack 2 (A, D, I)
    make_move(current, 3, 2, moves)  # Move B
    make_move(current, 3, 2, moves)  # Move D
    make_move(current, 3, 2, moves)  # Move C
    make_move(current, 3, 2, moves)  # Move I
    make_move(current, 3, 2, moves)  # Move A
    
    # Step 3: Build stack 3 (B, C, F)
    make_move(current, 2, 3, moves)  # Move I
    make_move(current, 2, 3, moves)  # Move C
    make_move(current, 2, 3, moves)  # Move B
    make_move(current, 1, 3, moves)  # Move F
    
    # Step 4: Complete stack 2
    make_move(current, 3, 2, moves)  # Move I

    return moves

# Initial and goal states
initial = [['F'], ['G', 'B', 'D'], ['A', 'I', 'C', 'E', 'H']]
goal = [['E', 'G', 'H'], ['A', 'D', 'I'], ['B', 'C', 'F']]

# Find and print solution
solution = solve_blocksworld(initial, goal)
print('<<<')
print('\n'.join(solution))
print('>>>')