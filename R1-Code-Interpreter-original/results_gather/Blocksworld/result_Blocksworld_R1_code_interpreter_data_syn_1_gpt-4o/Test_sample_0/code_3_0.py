# Initial configuration
stacks = {
    'stack1': ['E', 'B'],
    'stack2': ['A', 'D', 'C'],
    'stack3': []
}

# Revised moves to achieve the goal state
moves = [
    ('C', 'stack2', 'stack3'),
    ('D', 'stack2', 'stack1'),
    ('A', 'stack2', 'stack3'),
    ('D', 'stack1', 'stack3'),
    ('B', 'stack1', 'stack3'),
    ('E', 'stack1', 'stack3')
]

# Execute the moves
for block, source, destination in moves:
    stacks[destination].append(stacks[source].pop())

# Print the final configuration of stacks
print(stacks)