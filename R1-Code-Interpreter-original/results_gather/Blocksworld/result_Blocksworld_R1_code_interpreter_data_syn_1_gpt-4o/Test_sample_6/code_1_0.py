# Initial state
stack1 = ['A', 'C']
stack2 = ['E', 'D', 'B', 'H', 'F']
stack3 = ['G']

# Moves
moves = []

# Move G from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move G from 3 to 1")

# Move H from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move H from 2 to 3")

# Move B from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move B from 2 to 3")

# Move C from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move C from 1 to 3")

# Move D from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move D from 2 to 3")

# Move E from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move E from 2 to 3")

# Move H from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move H from 3 to 1")

# Print the moves
print("\n".join(moves))