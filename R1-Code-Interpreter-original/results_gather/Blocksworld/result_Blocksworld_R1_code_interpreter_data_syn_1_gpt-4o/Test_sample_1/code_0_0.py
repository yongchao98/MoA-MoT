# Initial and goal states
initial_state = {
    'stack1': [],
    'stack2': ['A', 'B'],
    'stack3': ['D'],
    'stack4': ['C', 'E']
}

goal_state = {
    'stack1': [],
    'stack2': ['A', 'D', 'E'],
    'stack3': ['B', 'C']
}

# Moves to achieve the goal state
moves = [
    "Move E from stack4 to stack1",
    "Move C from stack4 to stack3",
    "Move D from stack3 to stack2",
    "Move B from stack2 to stack3",
    "Move E from stack1 to stack2"
]

# Print the moves
for move in moves:
    print(move)