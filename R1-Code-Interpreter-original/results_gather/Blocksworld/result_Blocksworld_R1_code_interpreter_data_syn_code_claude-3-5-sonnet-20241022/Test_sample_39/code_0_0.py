def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")
    print()

def find_solution():
    # Initial state
    stacks = [
        ['H', 'B', 'F'],
        ['I'],
        ['G', 'C', 'D'],
        ['J', 'A', 'E']
    ]
    
    # Goal state
    goal = [
        ['A', 'C', 'D', 'F', 'G'],
        ['E', 'J'],
        ['B', 'H', 'I']
    ]
    
    # Manual solution (based on analysis)
    moves = [
        (1, 2),  # Move F from 1 to 2
        (1, 3),  # Move B from 1 to 3
        (1, 2),  # Move H from 1 to 2
        (4, 1),  # Move E from 4 to 1
        (4, 2),  # Move A from 4 to 2
        (4, 1),  # Move J from 4 to 1
        (3, 1),  # Move D from 3 to 1
        (3, 4),  # Move C from 3 to 4
        (3, 1),  # Move G from 3 to 1
        (2, 3),  # Move H from 2 to 3
        (2, 4),  # Move I from 2 to 4
        (2, 3),  # Move F from 2 to 3
        (2, 1),  # Move A from 2 to 1
        (4, 2),  # Move I from 4 to 2
        (4, 2),  # Move C from 4 to 2
        (3, 1),  # Move F from 3 to 1
        (3, 1),  # Move H from 3 to 1
        (2, 1),  # Move C from 2 to 1
        (2, 1),  # Move I from 2 to 1
        (1, 2),  # Move I from 1 to 2
        (1, 2),  # Move H from 1 to 2
        (1, 3),  # Move C from 1 to 3
        (2, 3),  # Move H from 2 to 3
        (2, 1),  # Move I from 2 to 1
        (1, 2),  # Move I from 1 to 2
        (1, 2),  # Move F from 1 to 2
        (2, 1),  # Move F from 2 to 1
        (2, 1),  # Move I from 2 to 1
    ]
    
    solution = []
    for from_stack, to_stack in moves:
        if is_valid_move(stacks, from_stack, to_stack):
            block = make_move(stacks, from_stack, to_stack)
            solution.append(f"Move {block} from {from_stack} to {to_stack}")
    
    print("\n".join(solution))

find_solution()