# Initial State
stack1 = ['A', 'E']
stack2 = ['D']
stack3 = ['C', 'B']

# Moves
moves = []

# Move B from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move B from 3 to 1")

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from 3 to 2")

# Move D from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move D from 2 to 3")

# Print the moves
for move in moves:
    print(move)