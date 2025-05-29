def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def get_top_blocks(stacks):
    return [stack[-1] if stack else None for stack in stacks]

def is_goal_reached(current, goal):
    return current == goal

def solve_blocksworld():
    # Initial state
    initial_state = [
        ['A', 'G'],
        ['F', 'E'],
        ['C'],
        ['H', 'J', 'I', 'B', 'D']
    ]
    
    # Goal state
    goal_state = [
        ['A', 'B'],
        ['C', 'E', 'G', 'J'],
        ['D', 'F', 'H', 'I']
    ]
    
    moves = []
    current_state = [stack[:] for stack in initial_state]
    
    # Strategy: First clear the blocks that are not in their final position
    # Move D from stack4 to stack3
    moves.append(make_move(current_state, 4, 3))
    
    # Move B from stack4 to stack1
    moves.append(make_move(current_state, 4, 1))
    
    # Move I from stack4 to stack3
    moves.append(make_move(current_state, 4, 3))
    
    # Move J from stack4 to stack2
    moves.append(make_move(current_state, 4, 2))
    
    # Move H from stack4 to stack3
    moves.append(make_move(current_state, 4, 3))
    
    # Move G from stack1 to stack2
    moves.append(make_move(current_state, 1, 2))
    
    # Move E from stack2 to stack2 (it's already in stack2)
    
    # Move F from stack2 to stack3
    moves.append(make_move(current_state, 2, 3))
    
    # Move C from stack3 to stack2
    moves.append(make_move(current_state, 3, 2))
    
    print('\n'.join(moves))

solve_blocksworld()