# Define the initial state
stacks = {
    1: ['B', 'K', 'J', 'C'],
    2: ['F', 'H'],
    3: ['D', 'E', 'G', 'I', 'A']
}

# Define the correct moves to reach the goal state
moves = [
    ('A', 3, 2),
    ('I', 3, 1),
    ('D', 3, 1),
    ('F', 2, 1),
    ('J', 1, 2),
    ('B', 1, 2)
]

# Execute the moves
for block, source, destination in moves:
    stacks[destination].append(stacks[source].pop())

# Print the final state of the stacks
print(stacks)