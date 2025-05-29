# Define the initial and goal states
initial_state = {
    "stack1": ["C", "F"],
    "stack2": ["E", "B", "G", "H"],
    "stack3": ["A", "D"]
}

goal_state = {
    "stack1": ["A", "B", "F"],
    "stack2": ["H"],
    "stack3": ["C", "D", "E", "G"]
}

# Define the moves to reach the goal state
moves = [
    "Move D from stack3 to stack2",
    "Move A from stack3 to stack1",
    "Move D from stack2 to stack3",
    "Move H from stack2 to stack3",
    "Move G from stack2 to stack3",
    "Move B from stack2 to stack1",
    "Move E from stack2 to stack3",
    "Move C from stack1 to stack3",
    "Move F from stack1 to stack2",
    "Move A from stack1 to stack2",
    "Move B from stack1 to stack2"
]

# Print the moves
for move in moves:
    print(move)