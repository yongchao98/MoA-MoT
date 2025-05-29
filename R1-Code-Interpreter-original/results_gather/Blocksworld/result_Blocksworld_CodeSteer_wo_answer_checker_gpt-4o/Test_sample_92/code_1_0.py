from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Calculate the heuristic based on the number of misplaced blocks and their distances
    misplaced = 0
    for i, stack in enumerate(state):
        for j, block in enumerate(stack):
            if j >= len(goal_state[i]) or block != goal_state[i][j]:
                misplaced += 1
    return misplaced

def get_possible_moves(state, last_move):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j and (last_move is None or last_move != (j, i)):  # Move to a different stack and avoid reversing the last move
                    new_state = [list(s) for s in state]  # Deep copy the state
                    new_state[i].pop()  # Remove the block from the current stack
                    new_state[j].append(block)  # Add the block to the new stack
                    moves.append((tuple(tuple(s) for s in new_state), f"Move {block} from {i+1} to {j+1}", (i, j)))
    return moves

def a_star(initial_state, goal_state):
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

        for new_state, move, new_last_move in get_possible_moves(current_state, last_move):
            if new_state not in visited:
                cost = len(path) + 1 + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, path + [move], new_last_move))

    return None

initial_state = (('F',), ('G', 'B', 'D'), ('A', 'I', 'C', 'E', 'H'))
goal_state = (('E', 'G', 'H'), ('A', 'D', 'I'), ('B', 'C', 'F'))

solution = a_star(initial_state, goal_state)

if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")