from heapq import heappop, heappush

def heuristic(state, goal_state):
    misplaced = 0
    for i in range(3):
        for j in range(len(state[i])):
            if j < len(goal_state[i]) and state[i][j] != goal_state[i][j]:
                misplaced += 1
            elif j >= len(goal_state[i]):
                misplaced += 1
    return misplaced

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