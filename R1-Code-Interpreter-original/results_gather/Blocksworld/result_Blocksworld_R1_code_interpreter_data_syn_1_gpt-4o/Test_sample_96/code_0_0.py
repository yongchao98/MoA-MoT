# Initial state
stack1 = ['E', 'B', 'J', 'I']
stack2 = ['A', 'F', 'K', 'C']
stack3 = ['D', 'G', 'H']

# Moves to achieve the goal state
moves = [
    ('I', 'stack1', 'stack3'),
    ('C', 'stack2', 'stack1'),
    ('K', 'stack2', 'stack1'),
    ('F', 'stack2', 'stack1'),
    ('A', 'stack2', 'stack1'),
    ('E', 'stack1', 'stack3'),
    ('B', 'stack1', 'stack3'),
    ('H', 'stack3', 'stack2'),
    ('G', 'stack3', 'stack2'),
    ('D', 'stack3', 'stack2')
]

# Function to perform the moves
def perform_moves(moves, stack1, stack2, stack3):
    stacks = {'stack1': stack1, 'stack2': stack2, 'stack3': stack3}
    for block, source, destination in moves:
        if stacks[source][-1] == block:
            stacks[destination].append(stacks[source].pop())
            print(f"Move {block} from {source} to {destination}")

perform_moves(moves, stack1, stack2, stack3)