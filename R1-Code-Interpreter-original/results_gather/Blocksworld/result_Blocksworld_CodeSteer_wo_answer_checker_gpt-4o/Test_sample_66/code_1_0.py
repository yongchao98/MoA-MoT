def heuristic(stacks, goal):
    # Improved heuristic: sum of distances of misplaced blocks to their target positions
    distance = 0
    for stack, goal_stack in zip(stacks, goal):
        for i, block in enumerate(stack):
            if i < len(goal_stack) and block != goal_stack[i]:
                # Calculate the distance to the target position
                target_index = goal_stack.index(block) if block in goal_stack else len(goal_stack)
                distance += abs(i - target_index)
    return distance

def ida_star(initial_state, goal_state, time_limit=10):
    start_time = time.time()

    def search(path, g, bound):
        current_state = path[-1]
        f = g + heuristic(current_state, goal_state)
        if f > bound:
            return f
        if is_goal_state(current_state, goal_state):
            return path
        min_bound = float('inf')
        for move in get_possible_moves(current_state):
            new_state = apply_move(current_state, move)
            if new_state not in path:  # Avoid cycles
                path.append(new_state)
                result = search(path, g + 1, bound)
                if isinstance(result, list):
                    return result
                if result < min_bound:
                    min_bound = result
                path.pop()
        return min_bound

    bound = heuristic(initial_state, goal_state)
    path = [initial_state]

    while True:
        if time.time() - start_time > time_limit:
            print("Timeout reached, no solution found within time limit.")
            return None
        result = search(path, 0, bound)
        if isinstance(result, list):
            return result
        if result == float('inf'):
            return None
        bound = result

initial_state = [['B', 'A'], ['C', 'F', 'G'], ['D', 'I', 'E', 'J', 'H']]
goal_state = [['A', 'E', 'I'], ['B', 'C', 'J'], ['D', 'F', 'G', 'H']]

solution = ida_star(initial_state, goal_state)

if solution:
    print("<<<")
    for i in range(1, len(solution)):
        prev_state = solution[i-1]
        current_state = solution[i]
        for j, (prev_stack, current_stack) in enumerate(zip(prev_state, current_state)):
            if len(prev_stack) > len(current_stack):
                block = prev_stack[-1]
                src = j
            elif len(prev_stack) < len(current_stack):
                dest = j
        print(f"Move {block} from stack{src+1} to stack{dest+1}")
    print(">>>")
else:
    print("No solution found.")