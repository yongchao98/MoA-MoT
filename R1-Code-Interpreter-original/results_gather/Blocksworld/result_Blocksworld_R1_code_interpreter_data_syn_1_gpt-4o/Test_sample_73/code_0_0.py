# Define the initial and goal states
initial_state = {
    "stack1": ["D"],
    "stack2": ["F", "A", "B"],
    "stack3": ["E"],
    "stack4": ["C"]
}

goal_state = {
    "stack1": [],
    "stack2": ["B", "C", "D", "E"],
    "stack3": ["A", "F"]
}

# Define the moves to reach the goal state
moves = [
    "Move D from stack1 to stack2",
    "Move C from stack4 to stack2",
    "Move E from stack3 to stack2",
    "Move A from stack2 to stack3",
    "Move F from stack2 to stack3"
]

# Print the moves
for move in moves:
    print(move)