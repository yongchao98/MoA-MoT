# Define the initial and goal states
initial_state = {
    'stack1': ['E', 'C', 'J'],
    'stack2': ['H', 'D'],
    'stack3': ['G', 'I', 'A', 'B', 'F']
}

goal_state = {
    'stack1': ['A', 'C', 'D', 'E', 'F', 'G'],
    'stack2': ['B', 'H', 'J'],
    'stack3': ['I']
}

# Define the moves to reach the goal state
moves = [
    "Move F from stack3 to stack1",
    "Move B from stack3 to stack2",
    "Move A from stack3 to stack1",
    "Move G from stack3 to stack1",
    "Move D from stack2 to stack1",
    "Move J from stack1 to stack2"
]

# Print the moves
for move in moves:
    print(move)