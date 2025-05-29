# Define the initial and goal states
initial_state = {
    'stack1': ['G'],
    'stack2': ['F', 'A', 'C'],
    'stack3': ['B', 'D', 'E']
}

goal_state = {
    'stack1': ['A', 'D', 'G'],
    'stack2': [],
    'stack3': ['B', 'C', 'E', 'F']
}

# Define the moves to reach the goal state
moves = [
    "Move C from stack2 to stack3",
    "Move A from stack2 to stack1",
    "Move F from stack2 to stack3",
    "Move E from stack3 to stack2",
    "Move D from stack3 to stack1",
    "Move E from stack2 to stack3"
]

# Print the moves
for move in moves:
    print(move)