# Initial state
stack1 = []
stack2 = ['B']
stack3 = ['C', 'A', 'E', 'D']

# Moves
moves = []

# Move D from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move D from 3 to 1")

# Move E from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move E from 3 to 2")

# Move A from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move A from 3 to 1")

# Move E from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move E from 2 to 3")

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from 3 to 2")

# Move B from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move B from 2 to 3")

# Print the moves
print("\n".join(moves))