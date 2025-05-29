# Initial state
stacks = {
    1: ['G', 'F'],
    2: ['E'],
    3: ['A', 'B'],
    4: ['D', 'C']
}

# Moves to achieve the goal state
moves = [
    ('F', 1, 3),  # Move F from stack1 to stack3
    ('G', 1, 2),  # Move G from stack1 to stack2
    ('E', 2, 1),  # Move E from stack2 to stack1
    ('C', 4, 3),  # Move C from stack4 to stack3
    ('D', 4, 1),  # Move D from stack4 to stack1
    ('B', 3, 1),  # Move B from stack3 to stack1
    ('F', 3, 3),  # Move F from stack3 to stack3
    ('A', 3, 3)   # Move A from stack3 to stack3
]

# Execute the moves
for block, source, destination in moves:
    stacks[destination].append(stacks[source].pop())

# Print the final stacks
print(stacks)