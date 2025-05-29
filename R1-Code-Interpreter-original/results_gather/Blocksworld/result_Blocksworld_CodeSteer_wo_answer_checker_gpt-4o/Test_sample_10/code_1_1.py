import heapq

def heuristic(state, goal):
    # Count the number of blocks not in the correct position
    misplaced = 0
    for stack, goal_stack in zip(state, goal):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def get_successors(state):
    successors = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()
                    new_state[j].append(block)
                    successors.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return successors

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    g_costs = {initial_state: 0}
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)

        if current_state == goal_state:
            return path

        visited.add(current_state)

        for successor, move in get_successors(current_state):
            if successor in visited:
                continue

            tentative_g_cost = g_costs[current_state] + 1

            if successor not in g_costs or tentative_g_cost < g_costs[successor]:
                g_costs[successor] = tentative_g_cost
                f_cost = tentative_g_cost + heuristic(successor, goal_state)
                heapq.heappush(open_set, (f_cost, successor, path + [move]))

    return None

initial_state = (('J',), ('D', 'H', 'F', 'I', 'B', 'A', 'G'), ('E', 'C'))
goal_state = (('G', 'I'), ('A', 'B', 'C', 'D', 'E', 'F', 'H'), ('J',))

solution = a_star_search(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")