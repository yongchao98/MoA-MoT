# Define the initial and goal states
initial_state = {
    "stack1": ["E", "D", "C"],
    "stack2": ["F", "B"],
    "stack3": ["A"],
    "stack4": []
}

goal_state = {
    "stack1": ["E"],
    "stack2": ["C", "D"],
    "stack3": ["A", "B", "F"]
}

# Define the moves to reach the goal state
moves = [
    "Move C from stack1 to stack2",
    "Move D from stack1 to stack2",
    "Move B from stack2 to stack3",
    "Move F from stack2 to stack3"
]

# Print the moves
for move in moves:
    print(move)