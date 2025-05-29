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
    
    # First, get E, G, H to stack1
    # Move H
    if get_top_block(current[2]) == 'H':
        make_move(current, 3, 1, moves)
    
    # Move E
    if get_top_block(current[2]) == 'E':
        make_move(current, 3, 1, moves)
    
    # Clear path for G and move it
    while get_top_block(current[1]) != 'G':
        make_move(current, 2, 3, moves)
    make_move(current, 2, 1, moves)
    
    # Get A, D, I to stack2
    # Move blocks to clear path for A
    while get_top_block(current[2]) != 'A':
        make_move(current, 3, 2, moves)
    make_move(current, 3, 2, moves)  # Move A
    
    # Move D
    if get_top_block(current[1]) == 'D':
        make_move(current, 2, 2, moves)
    elif get_top_block(current[2]) == 'D':
        make_move(current, 3, 2, moves)
    
    # Move I
    if get_top_block(current[2]) == 'I':
        make_move(current, 3, 2, moves)
    
    # Get B, C, F to stack3
    # Move B
    if get_top_block(current[1]) == 'B':
        make_move(current, 2, 3, moves)
    
    # Move C
    if get_top_block(current[2]) == 'C':
        make_move(current, 3, 3, moves)
    
    # Move F
    if get_top_block(current[0]) == 'F':
        make_move(current, 1, 3, moves)
    
    return moves

# Initial and goal states
initial = [['F'], ['G', 'B', 'D'], ['A', 'I', 'C', 'E', 'H']]
goal = [['E', 'G', 'H'], ['A', 'D', 'I'], ['B', 'C', 'F']]

# Find and print solution
solution = solve_blocksworld(initial, goal)
print('<<<')
print('\n'.join(solution))
print('>>>')