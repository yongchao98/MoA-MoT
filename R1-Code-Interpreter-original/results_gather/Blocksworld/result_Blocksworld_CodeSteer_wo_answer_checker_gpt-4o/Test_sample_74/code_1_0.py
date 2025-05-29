import heapq

def heuristic(state, goal_state):
    misplaced = 0
    for i, stack in enumerate(state):
        goal_stack = goal_state[i]
        for j, block in enumerate(stack):
            if j >= len(goal_stack) or block != goal_stack[j]:
                misplaced += 1
    return misplaced

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((tuple(map(tuple, new_state)), f"Move {block} from stack{i+1} to stack{j+1}"))
    return moves

def blocks_world_astar(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0 + heuristic(initial_state, goal_state), 0, initial_state, []))
    visited = set()

    while open_set:
        _, cost, current_state, path = heapq.heappop(open_set)
        if current_state in visited:
            continue
        if current_state == goal_state:
            return path
        visited.add(current_state)
        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                new_cost = cost + 1
                heapq.heappush(open_set, (new_cost + heuristic(new_state, goal_state), new_cost, new_state, path + [move]))

initial_state = (('B', 'G', 'F'), ('A', 'E', 'C'), ('D',))
goal_state = (('D', 'F', 'G'), ('E',), ('A', 'B', 'C'))

solution = blocks_world_astar(initial_state, goal_state)
print("<<<" + "\n".join(solution) + ">>>")