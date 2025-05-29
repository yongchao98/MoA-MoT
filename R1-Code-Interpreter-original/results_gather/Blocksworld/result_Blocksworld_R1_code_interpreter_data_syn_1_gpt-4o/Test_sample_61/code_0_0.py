# Define the initial and goal states
initial_state = {
    'stack1': [],
    'stack2': ['B'],
    'stack3': ['G', 'D', 'F', 'A', 'C'],
    'stack4': ['E']
}

goal_state = {
    'stack1': ['A', 'B', 'D', 'F'],
    'stack2': ['C', 'E'],
    'stack3': ['G']
}

# Define the moves
moves = [
    "Move C from stack3 to stack2",
    "Move A from stack3 to stack1",
    "Move F from stack3 to stack1",
    "Move D from stack3 to stack1",
    "Move B from stack2 to stack1",
    "Move E from stack4 to stack2"
]

# Print the moves
for move in moves:
    print(move)