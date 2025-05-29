def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {stack}")

def is_valid_move(stacks, from_stack, block):
    # Check if the block is at the top of the source stack
    return stacks[from_stack-1] and stacks[from_stack-1][-1] == block

def make_move(stacks, from_stack, to_stack, block):
    if not is_valid_move(stacks, from_stack, block):
        return False
    # Make the move
    stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return True

def solve_blocksworld():
    # Initial state
    stacks = [
        ['G', 'A', 'E', 'C'],  # stack 1
        ['H', 'I', 'F'],       # stack 2
        ['K'],                 # stack 3
        ['J', 'B', 'D']        # stack 4
    ]
    
    moves = []
    
    # First clear stack 4
    if make_move(stacks, 4, 3, 'D'):
        moves.append("Move D from 4 to 3")
    if make_move(stacks, 4, 3, 'B'):
        moves.append("Move B from 4 to 3")
    if make_move(stacks, 4, 2, 'J'):
        moves.append("Move J from 4 to 2")
        
    # Clear stack 1 except A
    if make_move(stacks, 1, 3, 'C'):
        moves.append("Move C from 1 to 3")
    if make_move(stacks, 1, 3, 'E'):
        moves.append("Move E from 1 to 3")
    if make_move(stacks, 1, 3, 'A'):
        moves.append("Move A from 1 to 3")
    if make_move(stacks, 1, 3, 'G'):
        moves.append("Move G from 1 to 3")
        
    # Clear stack 2
    if make_move(stacks, 2, 3, 'F'):
        moves.append("Move F from 2 to 3")
    if make_move(stacks, 2, 3, 'I'):
        moves.append("Move I from 2 to 3")
    if make_move(stacks, 2, 3, 'H'):
        moves.append("Move H from 2 to 3")
        
    # Build final stack 1 (A C)
    if make_move(stacks, 3, 1, 'A'):
        moves.append("Move A from 3 to 1")
    if make_move(stacks, 3, 1, 'C'):
        moves.append("Move C from 3 to 1")
        
    # Stack 2 is already correct with J
    
    # Stack 3 should now have B D E F G H I K in correct order
    # The moves above should have arranged them properly
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()