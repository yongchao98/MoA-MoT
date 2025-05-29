# Initial state
stacks = {
    1: ['C', 'E'],
    2: [],
    3: ['A'],
    4: ['D', 'B']
}

# Moves to achieve the goal state
moves = []

# Move B from stack4 to stack2
block = stacks[4].pop()
stacks[2].append(block)
moves.append(f"Move {block} from 4 to 2")

# Move D from stack4 to stack1
block = stacks[4].pop()
stacks[1].append(block)
moves.append(f"Move {block} from 4 to 1")

# Move C from stack1 to stack2
block = stacks[1].pop()
stacks[2].append(block)
moves.append(f"Move {block} from 1 to 2")

# Print the moves
for move in moves:
    print(move)