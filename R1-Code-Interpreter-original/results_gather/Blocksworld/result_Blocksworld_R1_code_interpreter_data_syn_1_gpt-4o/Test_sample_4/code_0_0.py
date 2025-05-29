# Define the initial and goal states
initial_state = {
    'stack1': ['F', 'A', 'G'],
    'stack2': ['E', 'B'],
    'stack3': ['C', 'D']
}

goal_state = {
    'stack1': ['A', 'B', 'D'],
    'stack2': ['G'],
    'stack3': ['C', 'E', 'F']
}

# Define the moves to reach the goal state
moves = [
    "Move G from stack1 to stack2",
    "Move F from stack1 to stack3",
    "Move E from stack2 to stack3",
    "Move B from stack2 to stack1",
    "Move D from stack3 to stack1"
]

# Print the moves
for move in moves:
    print(move)