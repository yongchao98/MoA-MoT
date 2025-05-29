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

def solve_blocksworld():
    # Initial state
    current_state = [
        ['H', 'I', 'B', 'E', 'A'],
        ['J'],
        ['G'],
        ['F', 'D', 'K', 'C']
    ]
    
    # Goal state
    goal_state = [
        ['H'],
        ['B', 'C', 'F', 'G', 'I'],
        ['A', 'D', 'E', 'J', 'K'],
        []
    ]
    
    moves = []
    # Manual solution (based on analysis)
    # Move C from top of stack 4
    moves.append(make_move(current_state, 4, 2))  # Move C from 4 to 3
    # Move K from top of stack 4
    moves.append(make_move(current_state, 4, 3))  # Move K from 4 to 3
    # Move D from top of stack 4
    moves.append(make_move(current_state, 4, 3))  # Move D from 4 to 3
    # Move F from top of stack 4
    moves.append(make_move(current_state, 4, 2))  # Move F from 4 to 2
    # Move A from top of stack 1
    moves.append(make_move(current_state, 1, 3))  # Move A from 1 to 3
    # Move E from top of stack 1
    moves.append(make_move(current_state, 1, 3))  # Move E from 1 to 3
    # Move B from top of stack 1
    moves.append(make_move(current_state, 1, 2))  # Move B from 1 to 2
    # Move I from top of stack 1
    moves.append(make_move(current_state, 1, 2))  # Move I from 1 to 2
    # Move G from stack 3
    moves.append(make_move(current_state, 3, 2))  # Move G from 3 to 2
    # Move J from stack 2
    moves.append(make_move(current_state, 2, 3))  # Move J from 2 to 3

    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()