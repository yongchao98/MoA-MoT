from heapq import heappop, heappush

def heuristic(state, goal):
    misplaced = 0
    for stack, goal_stack in zip(state, goal):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def get_valid_moves(state):
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
    heappush(open_set, (0, initial_state, []))
    visited = set()
    visited.add(initial_state)

    while open_set:
        _, current_state, path = heappop(open_set)

        if current_state == goal_state:
            return path

        for new_state, move in get_valid_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                cost = len(path) + 1
                estimated_cost = cost + heuristic(new_state, goal_state)
                heappush(open_set, (estimated_cost, new_state, path + [move]))

initial_state = (('I', 'A', 'G'), ('D', 'F', 'E'), ('H', 'B', 'C'))
goal_state = (('B', 'C', 'H'), ('F', 'G', 'I'), ('A', 'D', 'E'))

solution = a_star(initial_state, goal_state)

print("<<<" + "\n".join(solution) + ">>>")