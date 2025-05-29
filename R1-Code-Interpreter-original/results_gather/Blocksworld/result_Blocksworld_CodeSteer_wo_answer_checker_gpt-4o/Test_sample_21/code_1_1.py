def move_blocks(initial_state, goal_state, moves=[]):
    if initial_state == goal_state:
        return moves

    for i, stack in enumerate(initial_state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(initial_state)):
                if i != j:  # Move to a different stack
                    # Create a new state by moving the block
                    new_state = [list(s) for s in initial_state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    # Recurse with the new state
                    result = move_blocks(new_state, goal_state, moves + [f"Move {block} from stack{i+1} to stack{j+1}"])
                    if result:
                        return result
    return None

initial_state = [[], ['F'], ['G', 'C'], ['A', 'B', 'D', 'E', 'H']]
goal_state = [['B', 'D', 'E'], ['G', 'H'], ['A', 'C', 'F'], []]

solution = move_blocks(initial_state, goal_state)
if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")