from heapq import heappop, heappush

def is_goal_state(state, goal):
    return state == goal

def heuristic(state, goal):
    # Simple heuristic: count the number of blocks in the correct position
    score = 0
    for stack, goal_stack in zip(state, goal):
        for block, goal_block in zip(stack, goal_stack):
            if block == goal_block:
                score += 1
            else:
                break
    return -score  # Negative because heapq is a min-heap

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [s[:] for s in state]  # Deep copy of the state
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

        if is_goal_state(current_state, goal_state):
            return path

        for new_state, move in get_possible_moves(current_state):
            if tuple(tuple(stack) for stack in new_state) not in visited:
                new_path = path + [move]
                priority = len(new_path) + heuristic(new_state, goal_state)
                heappush(open_set, (priority, new_state, new_path))

    return None

initial_state = [['A', 'E', 'F'], ['B', 'C'], ['D', 'H', 'G'], ['I']]
goal_state = [['B', 'C', 'D', 'E'], ['A', 'F', 'G', 'H'], ['I'], []]

solution = a_star_search(initial_state, goal_state)

if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")