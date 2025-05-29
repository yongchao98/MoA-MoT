def heuristic(state, goal_state):
    misplaced = 0
    out_of_order = 0
    for stack, goal_stack in zip(state, goal_state):
        for i, block in enumerate(stack):
            if i < len(goal_stack):
                if block != goal_stack[i]:
                    misplaced += 1
                if block in goal_stack and goal_stack.index(block) != i:
                    out_of_order += 1
            else:
                misplaced += 1
    return misplaced + out_of_order

def ida_star_blocks_world(initial_state, goal_state):
    def search(path, g, threshold):
        current_state = path[-1]
        f = g + heuristic(current_state, goal_state)
        if f > threshold:
            return f
        if current_state == goal_state:
            return path
        min_threshold = float('inf')
        for i, stack in enumerate(current_state):
            if stack:
                block = stack[-1]
                for j in range(len(current_state)):
                    if i != j:
                        new_state = list(map(list, current_state))
                        new_state[i].pop()
                        new_state[j].append(block)
                        new_state_tuple = tuple(map(tuple, new_state))
                        if new_state_tuple not in path:
                            path.append(new_state_tuple)
                            result = search(path, g + 1, threshold)
                            if isinstance(result, list):
                                return result
                            if result < min_threshold:
                                min_threshold = result
                            path.pop()
        return min_threshold

    threshold = heuristic(initial_state, goal_state)
    path = [initial_state]
    while True:
        result = search(path, 0, threshold)
        if isinstance(result, list):
            return result
        if result == float('inf'):
            return None
        threshold = result

# Initial and goal states
initial_state = ((), ('E', 'B', 'I'), ('F', 'D'), ('C', 'H', 'G', 'A'))
goal_state = (('A', 'B', 'D', 'H'), ('I',), ('C', 'E', 'F', 'G'))

# Find the sequence of moves
solution_path = ida_star_blocks_world(initial_state, goal_state)

# Extract and print the moves
if solution_path:
    moves = []
    for i in range(1, len(solution_path)):
        prev_state = solution_path[i-1]
        current_state = solution_path[i]
        for j, (prev_stack, current_stack) in enumerate(zip(prev_state, current_state)):
            if len(prev_stack) < len(current_stack):
                block = current_stack[-1]
                source = prev_state.index(tuple(prev_stack + (block,)))
                destination = j
                moves.append(f"Move {block} from {source+1} to {destination+1}")
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")