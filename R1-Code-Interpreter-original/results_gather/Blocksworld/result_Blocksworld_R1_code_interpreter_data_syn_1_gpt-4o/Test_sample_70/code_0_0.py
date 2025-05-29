# Initial State
stack1 = ['A']
stack2 = ['E', 'C']
stack3 = ['D', 'B']

# Moves
moves = []

# Move B from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move B from 3 to 1")

# Move D from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move D from 3 to 1")

# Move A from stack1 to stack3
stack3.append(stack1.pop(0))
moves.append("Move A from 1 to 3")

# Move C from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move C from 2 to 3")

# Move E from stack2 to stack2 (no move needed, already in place)
# But for the sake of completion, we can mention it
moves.append("Move E from 2 to 2")

# Print the moves
print("\n".join(moves))