def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def is_goal_state(stacks, goal):
    return all(current == target for current, target in zip(stacks, goal))

def solve_blocksworld():
    # Initial state
    current_stacks = [
        ['G', 'A', 'E', 'C'],
        ['H', 'I', 'F'],
        ['K'],
        ['J', 'B', 'D']
    ]
    
    # Goal state
    goal_stacks = [
        ['A', 'C'],
        ['J'],
        ['B', 'D', 'E', 'F', 'G', 'H', 'I', 'K']
    ]
    
    moves = []
    
    # Strategy: First clear stack3 to use as temporary storage
    # Then build the goal stacks one by one
    
    # Move K from stack3 to stack2
    moves.append(make_move(current_stacks, 3, 2))
    
    # Clear stack1 except A and C
    moves.append(make_move(current_stacks, 1, 3))  # Move C
    moves.append(make_move(current_stacks, 1, 3))  # Move E
    moves.append(make_move(current_stacks, 1, 3))  # Move A
    moves.append(make_move(current_stacks, 1, 3))  # Move G
    
    # Clear stack2
    moves.append(make_move(current_stacks, 2, 1))  # Move K
    moves.append(make_move(current_stacks, 2, 1))  # Move F
    moves.append(make_move(current_stacks, 2, 1))  # Move I
    moves.append(make_move(current_stacks, 2, 1))  # Move H
    
    # Clear stack4
    moves.append(make_move(current_stacks, 4, 2))  # Move D
    moves.append(make_move(current_stacks, 4, 2))  # Move B
    moves.append(make_move(current_stacks, 4, 2))  # Move J
    
    # Build goal stacks
    # First build stack3 (B D E F G H I K)
    moves.append(make_move(current_stacks, 2, 3))  # Move B
    moves.append(make_move(current_stacks, 2, 3))  # Move D
    moves.append(make_move(current_stacks, 3, 3))  # Move E
    moves.append(make_move(current_stacks, 1, 3))  # Move F
    moves.append(make_move(current_stacks, 3, 3))  # Move G
    moves.append(make_move(current_stacks, 1, 3))  # Move H
    moves.append(make_move(current_stacks, 1, 3))  # Move I
    moves.append(make_move(current_stacks, 1, 3))  # Move K
    
    # Build stack2 (J)
    moves.append(make_move(current_stacks, 2, 2))  # Move J
    
    # Build stack1 (A C)
    moves.append(make_move(current_stacks, 3, 1))  # Move A
    moves.append(make_move(current_stacks, 3, 1))  # Move C
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()