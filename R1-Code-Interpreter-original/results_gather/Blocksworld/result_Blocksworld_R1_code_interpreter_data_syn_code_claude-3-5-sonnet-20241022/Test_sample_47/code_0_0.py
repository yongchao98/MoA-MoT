def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if is_valid_move(stacks, from_stack, to_stack):
        block = stacks[from_stack-1].pop()
        stacks[to_stack-1].append(block)
        return True
    return False

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def get_top_block(stack):
    return stack[-1] if stack else None

def find_solution():
    # Initial state
    initial = [['D', 'G', 'E', 'C'], ['A'], ['F'], ['B']]
    goal = [['D', 'G'], ['A', 'C'], ['B', 'E', 'F'], []]
    
    moves = []
    current = initial.copy()
    
    # Step 1: Move C from stack1 to stack2
    if make_move(current, 1, 2):
        moves.append("Move C from 1 to 2")
    
    # Step 2: Move E from stack1 to stack3
    if make_move(current, 1, 3):
        moves.append("Move E from 1 to 3")
    
    # Step 3: Move B from stack4 to stack3
    if make_move(current, 4, 3):
        moves.append("Move B from 4 to 3")
    
    # Step 4: Move F from stack3 to stack3 (it's already in correct position)
    if make_move(current, 3, 3):
        moves.append("Move F from 3 to 3")
    
    print("\n".join(moves))

find_solution()