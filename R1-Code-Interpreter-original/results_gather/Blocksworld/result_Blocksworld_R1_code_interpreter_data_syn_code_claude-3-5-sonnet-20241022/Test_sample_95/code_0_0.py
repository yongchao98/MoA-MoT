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
    return all(stack == goal[i] for i, stack in enumerate(stacks))

def solve_blocksworld():
    # Initial state
    current = [
        ['J', 'C'],           # stack 1
        [],                   # stack 2
        ['H', 'G', 'A', 'E', 'D'],  # stack 3
        ['F', 'B', 'I']      # stack 4
    ]
    
    # Goal state
    goal = [
        ['B', 'H'],          # stack 1
        ['A', 'D', 'E', 'F', 'G', 'I'],  # stack 2
        ['C', 'J']           # stack 3
    ]
    
    moves = []
    
    # Step 1: Clear space in stack 2 for most blocks
    # Move I from stack 4
    if 'I' in current[3]:
        moves.append(make_move(current, 4, 2))
    
    # Move C from stack 1
    if 'C' in current[0]:
        moves.append(make_move(current, 1, 3))
    
    # Move J from stack 1
    if 'J' in current[0]:
        moves.append(make_move(current, 1, 3))
    
    # Move blocks from stack 3 to stack 2 in correct order
    while current[2] and current[2][-1] in ['D', 'E', 'A', 'G']:
        moves.append(make_move(current, 3, 2))
    
    # Move blocks from stack 4
    while current[3]:
        block = current[3][-1]
        if block == 'B':
            moves.append(make_move(current, 4, 1))
        else:
            moves.append(make_move(current, 4, 2))
    
    # Move H to complete stack 1
    if 'H' in current[2]:
        moves.append(make_move(current, 3, 1))
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

solve_blocksworld()