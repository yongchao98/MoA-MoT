# Define the initial and goal states
initial_state = {
    "stack1": [],
    "stack2": ["K", "D", "I", "G", "B", "J"],
    "stack3": ["H", "C", "E", "A", "F"]
}

goal_state = {
    "stack1": ["A", "E", "F"],
    "stack2": ["B", "G", "K"],
    "stack3": ["C", "D", "H", "I", "J"]
}

# Define the moves
moves = [
    "Move F from stack3 to stack1",
    "Move A from stack3 to stack1",
    "Move E from stack3 to stack1",
    "Move J from stack2 to stack3",
    "Move B from stack2 to stack2",
    "Move G from stack2 to stack2",
    "Move I from stack2 to stack3",
    "Move D from stack2 to stack3",
    "Move C from stack3 to stack3",
    "Move H from stack3 to stack3",
    "Move K from stack2 to stack2"
]

# Print the moves
for move in moves:
    print(move)