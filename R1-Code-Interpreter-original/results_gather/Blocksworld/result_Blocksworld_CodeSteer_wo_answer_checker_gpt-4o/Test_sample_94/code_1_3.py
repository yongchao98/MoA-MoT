import heapq

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
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]  # Deep copy
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((new_state, f"Move {block} from {i+1} to {j+1}"))
    return moves

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)

        if current_state == goal_state:
            return path

        state_tuple = tuple(tuple(stack) for stack in current_state)
        if state_tuple in visited:
            continue
        visited.add(state_tuple)

        for new_state, move in get_possible_moves(current_state):
            new_path = path + [move]
            cost = len(new_path) + heuristic(new_state, goal_state)
            heapq.heappush(open_set, (cost, new_state, new_path))

    return None

initial_state = [['J', 'B'], ['F', 'I', 'E'], ['D', 'A', 'C', 'H', 'G']]
goal_state = [['B', 'C', 'G', 'H', 'I', 'J'], ['D', 'F'], ['A', 'E']]

solution = a_star_search(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")