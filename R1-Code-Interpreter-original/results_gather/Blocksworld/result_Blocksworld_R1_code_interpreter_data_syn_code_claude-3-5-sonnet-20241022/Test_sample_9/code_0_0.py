def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)

def print_move(block, from_stack, to_stack):
    return f"Move {block} from {from_stack} to {to_stack}"

def is_goal_reached(current, goal):
    return current == goal

def solve_blocksworld():
    # Initial state
    current_state = [
        [],                  # stack1
        ['E', 'B', 'I'],    # stack2
        ['F', 'D'],         # stack3
        ['C', 'H', 'G', 'A']# stack4
    ]
    
    # Goal state
    goal_state = [
        ['A', 'B', 'D', 'H'],   # stack1
        ['I'],                   # stack2
        ['C', 'E', 'F', 'G']    # stack3
    ]
    
    moves = []
    
    # Move A to stack1 (it needs to be at bottom)
    moves.append(print_move('A', 4, 1))
    make_move(current_state, 4, 1)
    
    # Move G to stack3
    moves.append(print_move('G', 4, 3))
    make_move(current_state, 4, 3)
    
    # Move H to stack1
    moves.append(print_move('H', 4, 1))
    make_move(current_state, 4, 1)
    
    # Move C to stack3
    moves.append(print_move('C', 4, 3))
    make_move(current_state, 4, 3)
    
    # Move I to stack2
    moves.append(print_move('I', 2, 2))
    make_move(current_state, 2, 2)
    
    # Move B to stack1
    moves.append(print_move('B', 2, 1))
    make_move(current_state, 2, 1)
    
    # Move E to stack3
    moves.append(print_move('E', 2, 3))
    make_move(current_state, 2, 3)
    
    # Move D to stack1
    moves.append(print_move('D', 3, 1))
    make_move(current_state, 3, 1)
    
    # Move F to stack3
    moves.append(print_move('F', 3, 3))
    make_move(current_state, 3, 3)
    
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

solve_blocksworld()