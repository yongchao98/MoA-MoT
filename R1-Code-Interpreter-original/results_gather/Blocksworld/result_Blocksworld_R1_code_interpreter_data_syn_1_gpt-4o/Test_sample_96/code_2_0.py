# Initial state
stack1 = ['E', 'B', 'J', 'I']
stack2 = ['A', 'F', 'K', 'C']
stack3 = ['D', 'G', 'H']

# Revised moves to achieve the goal state
revised_moves = [
    ('I', 'stack1', 'stack3'),
    ('C', 'stack2', 'stack1'),
    ('K', 'stack2', 'stack1'),
    ('F', 'stack2', 'stack1'),
    ('A', 'stack2', 'stack1'),
    ('J', 'stack1', 'stack3'),
    ('E', 'stack1', 'stack3'),
    ('B', 'stack1', 'stack3'),
    ('H', 'stack3', 'stack2'),
    ('G', 'stack3', 'stack2'),
    ('D', 'stack3', 'stack2')
]

# Function to perform the revised moves
def perform_revised_moves(moves, stack1, stack2, stack3):
    stacks = {'stack1': stack1, 'stack2': stack2, 'stack3': stack3}
    for block, source, destination in moves:
        if stacks[source] and stacks[source][-1] == block:
            stacks[destination].append(stacks[source].pop())
            print(f"Move {block} from {source} to {destination}")

perform_revised_moves(revised_moves, stack1, stack2, stack3)