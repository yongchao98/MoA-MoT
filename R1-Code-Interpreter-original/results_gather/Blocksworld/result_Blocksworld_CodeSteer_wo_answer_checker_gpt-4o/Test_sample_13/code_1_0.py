from heapq import heappop, heappush

def heuristic(state, goal):
    # Count misplaced blocks and blocks out of order
    misplaced = 0
    out_of_order = 0
    for stack, goal_stack in zip(state, goal):
        for i, (block, goal_block) in enumerate(zip(stack, goal_stack)):
            if block != goal_block:
                misplaced += 1
            if i > 0 and stack[i] < stack[i-1]:
                out_of_order += 1
    return misplaced + out_of_order

def get_possible_moves(state, last_move):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j and (last_move is None or last_move != (j, i)):  # Avoid reversing the last move
                    new_state = [list(s) for s in state]  # Convert tuples to lists for manipulation
                    new_state[i].pop()  # Remove the block from the current stack
                    new_state[j].append(block)  # Add the block to the new stack
                    moves.append((tuple(tuple(s) for s in new_state), f"Move {block} from {i+1} to {j+1}", (i, j)))
    return moves

def a_star_search(initial_state, goal_state):
    open_set = []
    heappush(open_set, (0, initial_state, [], None))
    visited = set()

    while open_set:
        _, current_state, path, last_move = heappop(open_set)
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if current_state == goal_state:
            return path

        for new_state, move, move_tuple in get_possible_moves(current_state, last_move):
            if new_state not in visited:
                new_path = path + [move]
                cost = len(new_path) + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, new_path, move_tuple))

    return None

initial_state = (('K', 'I', 'G'), ('A', 'F', 'B'), ('J', 'E'), ('C', 'H', 'D'))
goal_state = (('B', 'E', 'F', 'H'), ('A', 'D', 'I'), ('C', 'G', 'J', 'K'))

solution = a_star_search(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")