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
    
    # Step 1: Move H and E to stack 1 (they're on top of stack 3)
    make_move(current, 3, 1, moves)  # Move H
    make_move(current, 3, 1, moves)  # Move E
    
    # Step 2: Clear stack 2 to get G
    make_move(current, 2, 3, moves)  # Move D
    make_move(current, 2, 3, moves)  # Move B
    make_move(current, 2, 1, moves)  # Move G
    
    # Step 3: Move blocks to get A, D, I in stack 2
    make_move(current, 3, 2, moves)  # Move B
    make_move(current, 3, 2, moves)  # Move D
    make_move(current, 3, 2, moves)  # Move C
    make_move(current, 3, 2, moves)  # Move I
    make_move(current, 3, 2, moves)  # Move A
    
    # Step 4: Move blocks to get B, C, F in stack 3
    make_move(current, 2, 3, moves)  # Move I
    make_move(current, 2, 3, moves)  # Move C
    make_move(current, 2, 3, moves)  # Move B
    make_move(current, 1, 3, moves)  # Move F
    
    # Step 5: Final adjustments for stack 2
    make_move(current, 3, 2, moves)  # Move I
    make_move(current, 2, 1, moves)  # Move D to temp
    make_move(current, 2, 1, moves)  # Move A to temp
    make_move(current, 1, 2, moves)  # Move A back
    make_move(current, 1, 2, moves)  # Move D back
    make_move(current, 3, 2, moves)  # Move I back
    
    return moves

# Initial and goal states
initial = [['F'], ['G', 'B', 'D'], ['A', 'I', 'C', 'E', 'H']]
goal = [['E', 'G', 'H'], ['A', 'D', 'I'], ['B', 'C', 'F']]

# Find and print solution
solution = solve_blocksworld(initial, goal)
print('<<<')
print('\n'.join(solution))
print('>>>')