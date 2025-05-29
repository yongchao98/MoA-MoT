# Initial state
stacks = {
    1: ['F'],
    2: ['B', 'D'],
    3: ['C', 'E', 'A'],
    4: []
}

# Moves to achieve the goal state
moves = [
    ('A', 3, 4),  # Move A from stack3 to stack4
    ('E', 3, 1),  # Move E from stack3 to stack1
    ('C', 3, 2),  # Move C from stack3 to stack2
    ('D', 2, 1),  # Move D from stack2 to stack1
    ('F', 1, 3),  # Move F from stack1 to stack3
    ('A', 4, 3)   # Move A from stack4 to stack3
]

# Execute moves
for block, source, destination in moves:
    stacks[destination].append(stacks[source].pop())

# Print the final stacks to verify
print(stacks)