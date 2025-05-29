# Initial state
stack1 = ['D']
stack2 = ['A']
stack3 = ['B', 'C', 'E']

# Moves to achieve the goal state
moves = []

# Move E from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move E from 3 to 1")

# Move C from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move C from 3 to 1")

# Move B from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move B from 3 to 1")

# Move D from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move D from 1 to 3")

# Move C from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move C from 1 to 3")

# Move B from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move B from 1 to 3")

# Move E from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move E from 1 to 3")

# Move D from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move D from 3 to 1")

# Move C from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move C from 3 to 1")

# Move B from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move B from 3 to 1")

# Move E from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move E from 3 to 1")

# Print the moves
print("\n".join(moves))