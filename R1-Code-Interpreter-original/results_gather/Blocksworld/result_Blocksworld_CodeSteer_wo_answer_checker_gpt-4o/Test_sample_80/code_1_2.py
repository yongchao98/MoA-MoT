import heapq

def is_goal(state, goal):
    return state == goal

def heuristic(state, goal):
    # Count the number of blocks not in their goal position
    misplaced = 0
    for stack, goal_stack in zip(state, goal):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
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
                    moves.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return moves

def a_star(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)
        if is_goal(current_state, goal_state):
            return path
        visited.add(current_state)
        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                cost = len(path) + 1 + heuristic(new_state, goal_state)
                heapq.heappush(open_set, (cost, new_state, path + [move]))
    return None

initial_state = (('E', 'C', 'J'), ('H', 'D'), ('G', 'I', 'A', 'B', 'F'))
goal_state = (('A', 'C', 'D', 'E', 'F', 'G'), ('B', 'H', 'J'), ('I',))

solution = a_star(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")