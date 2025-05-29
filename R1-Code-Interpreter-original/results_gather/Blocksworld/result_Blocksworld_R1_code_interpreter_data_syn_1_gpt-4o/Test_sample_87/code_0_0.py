# Define the sequence of moves to achieve the goal state
moves = [
    "Move F from stack4 to stack2",
    "Move E from stack4 to stack2",
    "Move A from stack4 to stack2",
    "Move C from stack3 to stack2",
    "Move B from stack1 to stack2",
    "Move D from stack1 to stack2",
    "Move B from stack2 to stack1",
    "Move D from stack2 to stack1"
]

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")