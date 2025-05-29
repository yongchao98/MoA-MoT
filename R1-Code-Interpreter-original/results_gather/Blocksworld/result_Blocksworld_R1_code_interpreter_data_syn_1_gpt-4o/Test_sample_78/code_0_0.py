# Initial State
stack1 = ['F', 'I', 'H']
stack2 = ['G', 'A', 'E', 'C', 'B']
stack3 = ['D']

# Goal State
goal_stack1 = ['A', 'C']
goal_stack2 = ['E', 'F', 'G']
goal_stack3 = ['B', 'D', 'H', 'I']

# Moves to achieve the goal state
moves = []

# Move B from stack2 to stack3
moves.append("Move B from 2 to 3")
stack3.append(stack2.pop())

# Move C from stack2 to stack1
moves.append("Move C from 2 to 1")
stack1.append(stack2.pop())

# Move E from stack2 to stack1
moves.append("Move E from 2 to 1")
stack1.append(stack2.pop())

# Move A from stack2 to stack1
moves.append("Move A from 2 to 1")
stack1.append(stack2.pop())

# Move G from stack2 to stack3
moves.append("Move G from 2 to 3")
stack3.append(stack2.pop())

# Move F from stack1 to stack2
moves.append("Move F from 1 to 2")
stack2.append(stack1.pop())

# Move I from stack1 to stack3
moves.append("Move I from 1 to 3")
stack3.append(stack1.pop())

# Move H from stack1 to stack3
moves.append("Move H from 1 to 3")
stack3.append(stack1.pop())

# Move D from stack3 to stack2
moves.append("Move D from 3 to 2")
stack2.append(stack3.pop())

# Move E from stack1 to stack2
moves.append("Move E from 1 to 2")
stack2.append(stack1.pop())

# Move A from stack1 to stack2
moves.append("Move A from 1 to 2")
stack2.append(stack1.pop())

# Print the moves
for move in moves:
    print(move)