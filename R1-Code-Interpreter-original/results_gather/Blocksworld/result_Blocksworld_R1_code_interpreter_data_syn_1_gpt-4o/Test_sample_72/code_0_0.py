# Initial state
stack1 = ['F']
stack2 = ['B']
stack3 = ['D', 'A', 'E', 'C']

# Moves to reach the goal state
moves = []

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from 3 to 2")

# Move E from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move E from 3 to 2")

# Move A from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move A from 3 to 1")

# Move B from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move B from 2 to 3")

# Move D from stack3 to stack3 (no actual move needed, just for clarity)
# stack3.append(stack3.pop()) # D is already in place, no move needed

# Move A from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move A from 1 to 3")

# Print the moves
print("\n".join(moves))