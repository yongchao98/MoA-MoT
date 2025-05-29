# Define the initial and goal states
initial_state = {
    'stack1': ['D'],
    'stack2': ['A', 'J'],
    'stack3': ['H', 'F', 'C', 'E', 'I'],
    'stack4': ['B', 'G']
}

goal_state = {
    'stack1': ['F', 'I'],
    'stack2': ['A', 'B', 'C', 'D', 'G'],
    'stack3': ['E', 'H', 'J']
}

# Define the moves to reach the goal state
moves = [
    "Move I from stack3 to stack1",
    "Move E from stack3 to stack3",
    "Move C from stack3 to stack2",
    "Move F from stack3 to stack1",
    "Move H from stack3 to stack3",
    "Move J from stack2 to stack3",
    "Move G from stack4 to stack2",
    "Move B from stack4 to stack2",
    "Move D from stack1 to stack2"
]

# Print the moves
for move in moves:
    print(move)