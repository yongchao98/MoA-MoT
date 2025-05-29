from heapq import heappop, heappush

def heuristic(state, goal):
    misplaced = 0
    for i, stack in enumerate(state):
        for j, block in enumerate(stack):
            if j < len(goal[i]) and block != goal[i][j]:
                misplaced += 1
            elif j >= len(goal[i]):
                misplaced += 1
    return misplaced

def get_neighbors(state):
    neighbors = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, state))
                    new_state[i] = new_state[i][:-1]
                    new_state[j].append(block)
                    neighbors.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return neighbors

def a_star(initial, goal):
    open_set = []
    heappush(open_set, (0 + heuristic(initial, goal), 0, initial, []))
    closed_set = set()

    while open_set:
        _, cost, current, path = heappop(open_set)

        if current == goal:
            return path

        if current in closed_set:
            continue

        closed_set.add(current)

        for neighbor, action in get_neighbors(current):
            if neighbor not in closed_set:
                heappush(open_set, (cost + 1 + heuristic(neighbor, goal), cost + 1, neighbor, path + [action]))

    return None

initial_state = ((), ('B', 'A', 'C'), ('E', 'D'))
goal_state = (('B', 'E'), ('C',), ('A', 'D'))

solution = a_star(initial_state, goal_state)
print("<<<" + "\n".join(solution) + ">>>")