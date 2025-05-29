from heapq import heappop, heappush

def heuristic(state, goal_state):
    misplaced = 0
    for stack, goal_stack in zip(state, goal_state):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def a_star(initial_state, goal_state):
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heappop(open_set)

        if current_state in visited:
            continue
        visited.add(current_state)

        if current_state == goal_state:
            return path

        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                new_path = path + [move]
                cost = len(new_path) + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, new_path))

    return None

initial_state = (('D',), ('G', 'C', 'A', 'F'), ('B', 'I', 'E', 'H'))
goal_state = (('A', 'F', 'H'), ('C', 'G', 'I'), ('B', 'D', 'E'))

solution = a_star(initial_state, goal_state)
print("<<<" + "\n".join(solution) + ">>>")