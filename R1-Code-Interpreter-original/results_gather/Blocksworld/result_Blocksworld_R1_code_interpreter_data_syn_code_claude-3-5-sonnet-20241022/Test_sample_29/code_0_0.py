def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if is_valid_move(stacks, from_stack, to_stack):
        block = stacks[from_stack].pop()
        stacks[to_stack].append(block)
        return True
    return False

def get_state_string(stacks):
    return str(stacks)

def find_solution():
    # Initialize stacks
    stacks = [
        ['D'],        # stack1
        ['C'],        # stack2
        ['E', 'B', 'A'],  # stack3
        []           # stack4
    ]
    
    goal_stacks = [
        ['D', 'E'],   # stack1
        ['C'],        # stack2
        ['A', 'B'],   # stack3
        []           # stack4
    ]
    
    moves = []
    
    # Move A from stack3 to stack4 (temporary)
    make_move(stacks, 2, 3)
    moves.append(f"Move A from 3 to 4")
    
    # Move B from stack3 to stack4
    make_move(stacks, 2, 3)
    moves.append(f"Move B from 3 to 4")
    
    # Move E from stack3 to stack1
    make_move(stacks, 2, 0)
    moves.append(f"Move E from 3 to 1")
    
    # Move B from stack4 to stack3
    make_move(stacks, 3, 2)
    moves.append(f"Move B from 4 to 3")
    
    # Move A from stack4 to stack3
    make_move(stacks, 3, 2)
    moves.append(f"Move A from 4 to 3")
    
    print("\n".join(moves))

find_solution()