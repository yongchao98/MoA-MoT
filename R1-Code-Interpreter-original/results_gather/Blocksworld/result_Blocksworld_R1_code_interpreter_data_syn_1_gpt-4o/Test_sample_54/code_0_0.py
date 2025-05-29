# Initial state
stacks = {
    1: ['F', 'B', 'I', 'K', 'A'],
    2: ['H', 'D'],
    3: ['J', 'E', 'C', 'G']
}

# Moves to achieve the goal state
moves = [
    ('G', 3, 1),
    ('C', 3, 1),
    ('E', 3, 2),
    ('J', 3, 3),
    ('D', 2, 2),
    ('H', 2, 3),
    ('A', 1, 2),
    ('K', 1, 2),
    ('I', 1, 2),
    ('B', 1, 2),
    ('F', 1, 3)
]

# Execute the moves
for block, source, destination in moves:
    stacks[destination].append(stacks[source].pop())

# Print the final state of the stacks
print(stacks)