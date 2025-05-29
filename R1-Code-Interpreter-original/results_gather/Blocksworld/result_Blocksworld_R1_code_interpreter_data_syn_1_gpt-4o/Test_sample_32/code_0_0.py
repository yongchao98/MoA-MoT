# Initial state
stack1 = ['G', 'D', 'B', 'C']
stack2 = ['F', 'A', 'E']
stack3 = []

# Moves
moves = []

# Move C from stack1 to stack2
stack2.append(stack1.pop())
moves.append("Move C from stack1 to stack2")

# Move E from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move E from stack2 to stack3")

# Move A from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move A from stack2 to stack1")

# Move F from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move F from stack2 to stack3")

# Move E from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move E from stack3 to stack1")

# Move F from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move F from stack3 to stack2")

# Move C from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move C from stack2 to stack3")

# Move F from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move F from stack2 to stack3")

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from stack3 to stack2")

# Print the moves
print("\n".join(moves))