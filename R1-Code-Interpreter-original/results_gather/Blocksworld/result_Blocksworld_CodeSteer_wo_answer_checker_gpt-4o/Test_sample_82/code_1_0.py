from heapq import heappop, heappush

def heuristic(state, goal):
    # Improved heuristic: count misplaced blocks and consider their order
    misplaced = 0
    for i, stack in enumerate(state):
        for j, block in enumerate(stack):
            if i < len(goal) and j < len(goal[i]) and block != goal[i][j]:
                misplaced += 1
    return misplaced

def get_successors(state):
    successors = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move block to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()
                    new_state[j].append(block)
                    successors.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return successors

def a_star(initial, goal):
    open_set = []
    heappush(open_set, (0, initial, []))
    visited = set()

    while open_set:
        _, current, path = heappop(open_set)

        if current == goal:
            return path

        if current in visited:
            continue

        visited.add(current)

        for successor, move in get_successors(current):
            if successor not in visited:
                new_path = path + [move]
                cost = len(new_path) + heuristic(successor, goal)
                heappush(open_set, (cost, successor, new_path))

    return None

initial_state = (('F', 'A'), ('K', 'D', 'E', 'B', 'I', 'C'), ('H', 'J', 'G'))
goal_state = (('G',), ('B', 'D'), ('A', 'C', 'E', 'F', 'H', 'I', 'J', 'K'))

solution = a_star(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")