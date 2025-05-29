# Initial configuration
stack1 = ['F', 'A']
stack2 = ['K', 'D', 'E', 'B', 'I', 'C']
stack3 = ['H', 'J', 'G']

# Moves to achieve the goal state
moves = []

# Move C from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move C from 2 to 3")

# Move I from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move I from 2 to 3")

# Move B from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move B from 2 to 3")

# Move E from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move E from 2 to 3")

# Move D from stack2 to stack2
stack2.append(stack2.pop())
moves.append("Move D from 2 to 2")

# Move K from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move K from 2 to 3")

# Move A from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move A from 1 to 3")

# Move F from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move F from 1 to 3")

# Move G from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move G from 3 to 1")

# Move B from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move B from 3 to 2")

# Move D from stack2 to stack2
stack2.append(stack2.pop())
moves.append("Move D from 2 to 2")

# Print the moves
print("\n".join(moves))