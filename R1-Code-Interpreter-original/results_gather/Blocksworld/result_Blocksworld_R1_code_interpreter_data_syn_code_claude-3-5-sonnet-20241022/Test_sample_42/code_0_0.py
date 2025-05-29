def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack].pop()
    stacks[to_stack].append(block)
    return f"Move {block} from {from_stack + 1} to {to_stack + 1}"

def is_goal_state(stacks, goal):
    return all(stack == goal[i] for i, stack in enumerate(stacks))

def solve_blocksworld():
    # Initial state
    stacks = [['A'], [], ['B', 'C', 'D', 'E']]
    
    # Goal state
    goal = [['D'], ['B', 'E'], ['A', 'C']]
    
    moves = []
    
    # Simple heuristic solution
    # First clear the way to get required blocks
    # Move E from stack3 to stack2
    moves.append(make_move(stacks, 2, 1))  # Move E to stack1
    # Move D from stack3 to stack1
    moves.append(make_move(stacks, 2, 0))  # Move D to goal position
    # Move E from stack1 to stack2
    moves.append(make_move(stacks, 1, 1))  # Move E to stack2
    # Move C from stack3 to stack1
    moves.append(make_move(stacks, 2, 0))  # Move C temporarily
    # Move B from stack3 to stack2
    moves.append(make_move(stacks, 2, 1))  # Move B to goal position
    # Move C from stack1 to stack3
    moves.append(make_move(stacks, 0, 2))  # Move C to goal position
    # Move A from stack1 to stack3
    moves.append(make_move(stacks, 0, 2))  # Move A to goal position
    
    print("<<<" + "\n".join(moves) + ">>>")

solve_blocksworld()