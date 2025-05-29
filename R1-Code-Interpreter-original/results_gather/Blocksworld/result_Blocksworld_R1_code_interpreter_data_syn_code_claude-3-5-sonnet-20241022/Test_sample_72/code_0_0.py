def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    # Initial state
    current_stacks = [['F'], ['B'], ['D', 'A', 'E', 'C']]
    goal_stacks = [['F'], ['C', 'E'], ['A', 'B', 'D']]
    
    # Store states and moves
    moves = []
    
    # Manual solution (based on analysis)
    # Move C from stack3 to stack2
    moves.append(make_move(current_stacks, 3, 2))
    
    # Move E from stack3 to stack1
    moves.append(make_move(current_stacks, 3, 1))
    
    # Move A from stack3 to stack1
    moves.append(make_move(current_stacks, 3, 1))
    
    # Move B from stack2 to stack3
    moves.append(make_move(current_stacks, 2, 3))
    
    # Move A from stack1 to stack3
    moves.append(make_move(current_stacks, 1, 3))
    
    # Move E from stack1 to stack2
    moves.append(make_move(current_stacks, 1, 2))
    
    print("<<<" + "\n".join(moves) + ">>>")

find_solution()