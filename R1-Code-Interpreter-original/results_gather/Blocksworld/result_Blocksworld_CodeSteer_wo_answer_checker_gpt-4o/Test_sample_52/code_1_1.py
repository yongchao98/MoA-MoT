def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source + 1} to {destination + 1}"

def is_goal_state(stacks, goal):
    return stacks == goal

def solve_blocksworld(initial, goal):
    stacks = [list(stack) for stack in initial]
    moves = []

    # Move A to stack2
    moves.append(move_block(stacks, 2, 1))
    # Move D to stack2
    moves.append(move_block(stacks, 2, 1))
    # Move H to stack1
    moves.append(move_block(stacks, 2, 0))
    # Move C to stack1
    moves.append(move_block(stacks, 1, 0))
    # Move E to stack3
    moves.append(move_block(stacks, 1, 2))
    # Move G to stack3
    moves.append(move_block(stacks, 1, 2))
    # Move B to stack1
    moves.append(move_block(stacks, 1, 0))
    # Move F to stack3
    moves.append(move_block(stacks, 0, 2))

    return moves

initial_state = [
    ['I', 'J', 'F'],
    ['B', 'G', 'E', 'C'],
    ['H', 'D', 'A']
]

goal_state = [
    ['B', 'C', 'H', 'I', 'J'],
    ['D'],
    ['A', 'E', 'F', 'G']
]

moves = solve_blocksworld(initial_state, goal_state)
print("<<<" + "\n".join(moves) + ">>>")