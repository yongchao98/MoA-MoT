from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Count the number of misplaced blocks compared to the goal state
    misplaced = 0
    for stack, goal_stack in zip(state, goal_state):
        for block, goal_block in zip(stack, goal_stack):
            if block != goal_block:
                misplaced += 1
    return misplaced

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]  # Deep copy of the state
                    new_state[i].pop()  # Remove the block from the current stack
                    new_state[j].append(block)  # Add the block to the new stack
                    moves.append((new_state, f"Move {block} from {i+1} to {j+1}"))
    return moves

def a_star_search(initial_state, goal_state):
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heappop(open_set)
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if current_state == goal_state:
            return path

        for new_state, move in get_possible_moves(current_state):
            if tuple(tuple(stack) for stack in new_state) not in visited:
                cost = len(path) + 1 + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, path + [move]))

    return None

initial_state = [['A', 'C'], ['E', 'D', 'B', 'H', 'F'], ['G']]
goal_state = [['A', 'G', 'H'], ['F'], ['B', 'C', 'D', 'E']]

solution = a_star_search(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")