def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, src, dest = move
    stacks[src].pop()
    stacks[dest].append(block)

def undo_move(stacks, move):
    block, src, dest = move
    stacks[dest].pop()
    stacks[src].append(block)

def solve(stacks, goal, moves=[]):
    if is_goal_state(stacks, goal):
        return moves

    for move in get_possible_moves(stacks):
        apply_move(stacks, move)
        result = solve(stacks, goal, moves + [move])
        if result:
            return result
        undo_move(stacks, move)

    return None

# Initial and goal states
initial_stacks = [['F'], ['B', 'D'], ['C', 'E', 'A'], []]
goal_stacks = [['D', 'E'], ['B', 'C'], ['A', 'F'], []]

# Solve the problem
solution = solve(initial_stacks, goal_stacks)

# Print the solution
if solution:
    print("<<<")
    for block, src, dest in solution:
        print(f"Move {block} from {src + 1} to {dest + 1}")
    print(">>>")
else:
    print("No solution found")