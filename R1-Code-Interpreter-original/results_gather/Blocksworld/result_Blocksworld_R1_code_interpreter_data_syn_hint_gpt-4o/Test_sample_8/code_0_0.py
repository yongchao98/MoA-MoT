# Initial state
stack1 = ['D']
stack2 = ['G', 'C', 'A', 'F']
stack3 = ['B', 'I', 'E', 'H']

# Moves to reach the goal state
moves = []

# Move H from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move H from 3 to 1")

# Move E from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move E from 3 to 2")

# Move I from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move I from 3 to 2")

# Move B from stack3 to stack3 (no actual move needed, just for clarity)
moves.append("Move B from 3 to 3")

# Move D from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move D from 1 to 3")

# Move F from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move F from 2 to 1")

# Move A from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move A from 2 to 1")

# Move G from stack2 to stack2 (no actual move needed, just for clarity)
moves.append("Move G from 2 to 2")

# Move C from stack2 to stack2 (no actual move needed, just for clarity)
moves.append("Move C from 2 to 2")

# Print the moves
print("\n".join(moves))