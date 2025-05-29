# Define the initial and goal states
initial_state = {
    "stack1": ["D", "C", "E"],
    "stack2": ["A"],
    "stack3": [],
    "stack4": ["B"]
}

goal_state = {
    "stack1": ["A", "C", "D"],
    "stack2": ["B"],
    "stack3": ["E"]
}

# Define the series of moves to reach the goal state
moves = [
    "Move E from stack1 to stack3",
    "Move D from stack1 to stack4",
    "Move A from stack2 to stack1",
    "Move D from stack4 to stack1",
    "Move B from stack4 to stack2"
]

# Print the moves
for move in moves:
    print(move)