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

def find_solution():
    # Initial state
    current_stacks = [
        ['F', 'H'],
        ['G', 'D'],
        ['E', 'C', 'B'],
        ['A', 'I']
    ]
    
    # Goal state
    goal_stacks = [
        ['D', 'F'],
        ['B', 'C', 'E'],
        ['A', 'G', 'H', 'I']
    ]
    
    moves = []
    # Manual solution (based on analysis)
    
    # First clear space and move blocks to achieve bottom layers
    moves.append(make_move(current_stacks, 1, 4))  # Move H from 1 to 4
    moves.append(make_move(current_stacks, 2, 1))  # Move D from 2 to 1
    moves.append(make_move(current_stacks, 3, 2))  # Move B from 3 to 2
    moves.append(make_move(current_stacks, 3, 2))  # Move C from 3 to 2
    moves.append(make_move(current_stacks, 3, 2))  # Move E from 3 to 2
    moves.append(make_move(current_stacks, 4, 3))  # Move I from 4 to 3
    moves.append(make_move(current_stacks, 4, 3))  # Move H from 4 to 3
    moves.append(make_move(current_stacks, 2, 3))  # Move G from 2 to 3
    moves.append(make_move(current_stacks, 4, 3))  # Move A from 4 to 3
    moves.append(make_move(current_stacks, 1, 2))  # Move F from 1 to 2
    moves.append(make_move(current_stacks, 1, 2))  # Move D from 1 to 2
    moves.append(make_move(current_stacks, 2, 1))  # Move D from 2 to 1
    moves.append(make_move(current_stacks, 2, 1))  # Move F from 2 to 1
    moves.append(make_move(current_stacks, 2, 3))  # Move E from 2 to 3
    moves.append(make_move(current_stacks, 2, 3))  # Move C from 2 to 3
    moves.append(make_move(current_stacks, 2, 3))  # Move B from 2 to 3
    moves.append(make_move(current_stacks, 3, 2))  # Move B from 3 to 2
    moves.append(make_move(current_stacks, 3, 2))  # Move C from 3 to 2
    moves.append(make_move(current_stacks, 3, 2))  # Move E from 3 to 2
    
    print("\n".join(moves))

find_solution()