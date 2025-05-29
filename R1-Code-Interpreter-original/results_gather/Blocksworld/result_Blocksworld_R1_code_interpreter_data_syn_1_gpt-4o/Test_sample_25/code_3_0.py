# Initial state
stacks = {
    1: ['A', 'G'],
    2: ['F', 'E'],
    3: ['C'],
    4: ['H', 'J', 'I', 'B', 'D']
}

# Corrected moves to achieve the goal state
moves = [
    ('D', 4, 3),
    ('B', 4, 1),
    ('I', 4, 3),
    ('J', 4, 2),
    ('H', 4, 3),
    ('G', 1, 2),
    ('F', 2, 3),
    ('C', 3, 2)
]

# Execute the moves
for block, source, destination in moves:
    if stacks[source][-1] == block:
        stacks[source].pop()
        stacks[destination].append(block)
    else:
        raise ValueError(f"Invalid move: {block} is not on top of stack {source}")

# Print the moves
for block, source, destination in moves:
    print(f"Move {block} from {source} to {destination}")