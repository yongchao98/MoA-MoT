# Define the sequence of moves
moves = [
    "Move C from stack1 to stack2",
    "Move B from stack1 to stack3",
    "Move A from stack1 to stack2",
    "Move F from stack3 to stack2",
    "Move A from stack2 to stack1",
    "Move C from stack2 to stack1"
]

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")