import heapq

def is_goal_state(state, goal_state):
    return state == goal_state

def get_valid_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(state, move):
    block, src, dest = move
    new_state = [stack[:] for stack in state]  # Deep copy the state
    new_state[src].pop()  # Remove the block from the source stack
    new_state[dest].append(block)  # Add the block to the destination stack
    return new_state

def state_to_tuple(state):
    return tuple(tuple(stack) for stack in state)

def heuristic(state, goal_state):
    misplaced_blocks = 0
    for stack, goal_stack in zip(state, goal_state):
        for i, block in enumerate(stack):
            if i < len(goal_stack) and block != goal_stack[i]:
                misplaced_blocks += 1
            elif i >= len(goal_stack):
                misplaced_blocks += 1
    return misplaced_blocks

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))  # (priority, state, path)
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)
        state_tuple = state_to_tuple(current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if is_goal_state(current_state, goal_state):
            return path

        for move in get_valid_moves(current_state):
            new_state = apply_move(current_state, move)
            new_path = path + [move]
            cost = len(new_path) + heuristic(new_state, goal_state)
            heapq.heappush(open_set, (cost, new_state, new_path))

    return None

initial_state = [['I', 'D', 'B'], ['G', 'J', 'A', 'C', 'E'], ['H', 'F']]
goal_state = [['E', 'H'], ['A', 'G'], ['B', 'C', 'D', 'F', 'I', 'J']]

solution = a_star_search(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(f"Move {block} from stack{src+1} to stack{dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")