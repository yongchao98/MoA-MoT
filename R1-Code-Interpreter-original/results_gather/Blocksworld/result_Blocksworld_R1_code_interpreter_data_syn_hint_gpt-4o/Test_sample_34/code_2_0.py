# Define the initial and goal states
initial_state = {
    'stack1': ['G', 'E', 'B'],
    'stack2': ['A'],
    'stack3': ['C', 'D', 'H', 'F']
}

goal_state = {
    'stack1': ['D', 'F', 'G'],
    'stack2': ['A', 'C', 'E', 'H'],
    'stack3': ['B']
}

# Define the moves to reach the goal state
moves = [
    "Move F from stack3 to stack1",
    "Move H from stack3 to stack2",
    "Move D from stack3 to stack1",
    "Move C from stack3 to stack2",
    "Move B from stack1 to stack3",
    "Move E from stack1 to stack2",
    "Move G from stack1 to stack1"
]

# Print the moves
for move in moves:
    print(move)