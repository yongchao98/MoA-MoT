from heapq import heappop, heappush

def heuristic(state, goal_state):
    distance = 0
    goal_positions = {}
    for i, stack in enumerate(goal_state):
        for j, block in enumerate(stack):
            goal_positions[block] = (i, j)

    for i, stack in enumerate(state):
        for j, block in enumerate(stack):
            if block in goal_positions:
                goal_i, goal_j = goal_positions[block]
                distance += abs(i - goal_i) + abs(j - goal_j)
    return distance

def a_star(initial_state, goal_state):
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = set()
    visited.add(tuple(map(tuple, initial_state)))

    while open_set:
        _, current_state, path = heappop(open_set)

        if current_state == goal_state:
            return path

        for i in range(3):
            if current_state[i]:  # If the stack is not empty
                block = current_state[i][-1]
                for j in range(3):
                    if i != j:  # Move to a different stack
                        new_state = [list(stack) for stack in current_state]
                        new_state[i].pop()
                        new_state[j].append(block)

                        state_tuple = tuple(map(tuple, new_state))
                        if state_tuple not in visited:
                            visited.add(state_tuple)
                            new_path = path + [f"Move {block} from {i+1} to {j+1}"]
                            cost = len(new_path) + heuristic(new_state, goal_state)
                            heappush(open_set, (cost, new_state, new_path))

    return None

initial_state = [['F', 'B', 'I', 'K', 'A'], ['H', 'D'], ['J', 'E', 'C', 'G']]
goal_state = [['C', 'G'], ['A', 'B', 'D', 'E', 'I', 'K'], ['F', 'H', 'J']]

solution = a_star(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")