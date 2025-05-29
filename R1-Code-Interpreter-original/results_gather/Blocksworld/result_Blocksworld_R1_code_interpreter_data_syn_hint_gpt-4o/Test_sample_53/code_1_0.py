# Define the initial and goal states
initial_state = {
    'stack1': ['C'],
    'stack2': ['B', 'E', 'A', 'J', 'F'],
    'stack3': ['G', 'I'],
    'stack4': ['D', 'H']
}

goal_state = {
    'stack1': ['E'],
    'stack2': ['A', 'B', 'C', 'D', 'H'],
    'stack3': ['F', 'G', 'I', 'J']
}

# Define the moves to achieve the goal state
moves = [
    "Move H from stack4 to stack2",
    "Move D from stack4 to stack2",
    "Move C from stack1 to stack2",
    "Move F from stack2 to stack3",
    "Move J from stack2 to stack3",
    "Move E from stack2 to stack1"
]

# Print the moves
for move in moves:
    print(move)