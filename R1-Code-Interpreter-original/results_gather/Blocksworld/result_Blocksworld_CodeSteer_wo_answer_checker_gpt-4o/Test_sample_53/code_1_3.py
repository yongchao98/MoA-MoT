import heapq

def is_goal_state(state, goal):
    return state == goal

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
    new_state = [list(stack) for stack in state]  # Deep copy the state
    new_state[src].pop()  # Remove block from source stack
    new_state[dest].append(block)  # Add block to destination stack
    return tuple(tuple(stack) for stack in new_state)

def heuristic(state, goal):
    # Calculate the sum of distances each block is from its goal position
    distance = 0
    goal_positions = {block: (i, j) for i, stack in enumerate(goal) for j, block in enumerate(stack)}
    for i, stack in enumerate(state):
        for j, block in enumerate(stack):
            if block in goal_positions:
                goal_i, goal_j = goal_positions[block]
                distance += abs(i - goal_i) + abs(j - goal_j)
    return distance

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = {}

    while open_set:
        _, current_state, moves = heapq.heappop(open_set)
        state_tuple = tuple(tuple(stack) for stack in current_state)
        
        if state_tuple in visited and visited[state_tuple] <= len(moves):
            continue
        visited[state_tuple] = len(moves)

        if is_goal_state(current_state, goal_state):
            return moves

        for move in get_valid_moves(current_state):
            new_state = apply_move(current_state, move)
            if new_state not in visited or visited[new_state] > len(moves) + 1:
                new_moves = moves + [move]
                cost = len(new_moves) + heuristic(new_state, goal_state)
                heapq.heappush(open_set, (cost, new_state, new_moves))

    return None

initial_state = (('C',), ('B', 'E', 'A', 'J', 'F'), ('G', 'I'), ('D', 'H'))
goal_state = (('E',), ('A', 'B', 'C', 'D', 'H'), ('F', 'G', 'I', 'J'))

solution_moves = a_star_search(initial_state, goal_state)

if solution_moves:
    print("<<<")
    for block, src, dest in solution_moves:
        print(f"Move {block} from {src + 1} to {dest + 1}")
    print(">>>")
else:
    print("No solution found.")