# Define the initial and goal states
initial_state = {
    "stack1": [],
    "stack2": ["F"],
    "stack3": ["G", "C"],
    "stack4": ["A", "B", "D", "E", "H"]
}

goal_state = {
    "stack1": ["B", "D", "E"],
    "stack2": ["G", "H"],
    "stack3": ["A", "C", "F"]
}

# Define the moves to reach the goal state
moves = [
    "Move H from stack4 to stack2",
    "Move E from stack4 to stack1",
    "Move D from stack4 to stack1",
    "Move B from stack4 to stack1",
    "Move A from stack4 to stack3",
    "Move C from stack3 to stack3",
    "Move G from stack3 to stack2",
    "Move F from stack2 to stack3"
]

# Print the moves
for move in moves:
    print(move)