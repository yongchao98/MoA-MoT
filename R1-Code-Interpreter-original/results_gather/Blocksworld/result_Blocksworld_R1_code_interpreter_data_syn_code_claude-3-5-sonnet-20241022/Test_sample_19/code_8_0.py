def print_state(state):
    for i, stack in enumerate(state):
        print(f"Stack {i+1}: {stack}")

def is_valid_move(state, block, from_stack, to_stack):
    # Check if block is at the top of from_stack
    if not state[from_stack] or state[from_stack][-1] != block:
        return False
    return True

def make_move(state, block, from_stack, to_stack):
    # Make a move and return new state
    new_state = [stack[:] for stack in state]
    new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

# Initial and goal states
initial = [['G', 'F'], ['E'], ['A', 'B'], ['D', 'C']]
goal = [['B', 'D', 'E'], ['G'], ['A', 'C', 'F']]

# Solution moves
moves = [
    ('F', 1, 4),      # Move F from stack1 to stack4
    ('G', 1, 3),      # Move G from stack1 to stack3 temporarily
    ('B', 3, 1),      # Move B from stack3 to stack1
    ('G', 3, 2),      # Move G to its final position in stack2
    ('F', 4, 3),      # Move F from stack4 to stack3
    ('C', 4, 3),      # Move C from stack4 to stack3
    ('D', 4, 1),      # Move D from stack4 to stack1
    ('E', 2, 1),      # Move E from stack2 to stack1
]

# Verify solution
current_state = [stack[:] for stack in initial]
solution = []
valid = True

for block, from_stack, to_stack in moves:
    if is_valid_move(current_state, block, from_stack-1, to_stack-1):
        current_state = make_move(current_state, block, from_stack-1, to_stack-1)
        solution.append(f"Move {block} from {from_stack} to {to_stack}")
    else:
        valid = False
        print(f"Invalid move: Move {block} from {from_stack} to {to_stack}")
        break

if valid:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("Solution is invalid")