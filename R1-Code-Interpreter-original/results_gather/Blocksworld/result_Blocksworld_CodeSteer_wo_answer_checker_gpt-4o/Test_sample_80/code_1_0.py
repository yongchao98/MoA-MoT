import heapq

def is_goal(state, goal):
    return state == goal

def heuristic(state, goal):
    misplaced = 0
    out_of_order = 0
    for stack, goal_stack in zip(state, goal):
        for i, (block, goal_block) in enumerate(zip(stack, goal_stack)):
            if block != goal_block:
                misplaced += 1
            if i > 0 and stack[i-1] != goal_stack[i-1]:
                out_of_order += 1
    return misplaced + out_of_order

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
                    moves.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return moves

def a_star(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = {}
    visited[initial_state] = 0

    while open_set:
        _, current_state, path = heapq.heappop(open_set)
        if is_goal(current_state, goal_state):
            return path
        for new_state, move in get_possible_moves(current_state):
            new_cost = len(path) + 1
            if new_state not in visited or new_cost < visited[new_state]:
                visited[new_state] = new_cost
                cost = new_cost + heuristic(new_state, goal_state)
                heapq.heappush(open_set, (cost, new_state, path + [move]))
    return None

initial_state = (('E', 'C', 'J'), ('H', 'D'), ('G', 'I', 'A', 'B', 'F'))
goal_state = (('A', 'C', 'D', 'E', 'F', 'G'), ('B', 'H', 'J'), ('I',))

solution = a_star(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")