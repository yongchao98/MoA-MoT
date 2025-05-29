def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def is_goal_state(stacks):
    goal = [['E'], ['A', 'B', 'C', 'F'], ['D'], []]
    return stacks == goal

def get_top_blocks(stacks):
    tops = []
    for stack in stacks:
        if stack:
            tops.append(stack[-1])
        else:
            tops.append(None)
    return tops

def find_solution():
    # Initial state
    stacks = [['A', 'D'], [], ['F', 'B'], ['E', 'C']]
    moves = []
    
    # Manual solution (based on analysis)
    # Move C from stack4 to stack2
    if is_valid_move(stacks, 4, 2):
        block = make_move(stacks, 4, 2)
        moves.append(f"Move {block} from 4 to 2")
    
    # Move B from stack3 to stack2
    if is_valid_move(stacks, 3, 2):
        block = make_move(stacks, 3, 2)
        moves.append(f"Move {block} from 3 to 2")
    
    # Move D from stack1 to stack3
    if is_valid_move(stacks, 1, 3):
        block = make_move(stacks, 1, 3)
        moves.append(f"Move {block} from 1 to 3")
    
    # Move E from stack4 to stack1
    if is_valid_move(stacks, 4, 1):
        block = make_move(stacks, 4, 1)
        moves.append(f"Move {block} from 4 to 1")
    
    # Move F from stack3 to stack2
    if is_valid_move(stacks, 3, 2):
        block = make_move(stacks, 3, 2)
        moves.append(f"Move {block} from 3 to 2")
    
    # Move A from stack1 to stack2
    if is_valid_move(stacks, 1, 2):
        block = make_move(stacks, 1, 2)
        moves.append(f"Move {block} from 1 to 2")
    
    print("<<<" + "\n".join(moves) + ">>>")

find_solution()