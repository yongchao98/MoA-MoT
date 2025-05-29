import heapq

def heuristic(state, goal):
    # Calculate the heuristic based on misplaced blocks and their distances
    misplaced = 0
    for stack, goal_stack in zip(state, goal):
        for i, block in enumerate(stack):
            if i >= len(goal_stack) or block != goal_stack[i]:
                misplaced += 1
    return misplaced

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]  # Deep copy
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((tuple(tuple(s) for s in new_state), f"Move {block} from {i+1} to {j+1}"))
    return moves

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)

        if current_state == goal_state:
            return path

        if current_state in visited:
            continue
        visited.add(current_state)

        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                new_path = path + [move]
                cost = len(new_path) + heuristic(new_state, goal_state)
                heapq.heappush(open_set, (cost, new_state, new_path))

    return None

initial_state = (('J', 'B'), ('F', 'I', 'E'), ('D', 'A', 'C', 'H', 'G'))
goal_state = (('B', 'C', 'G', 'H', 'I', 'J'), ('D', 'F'), ('A', 'E'))

solution = a_star_search(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")