from heapq import heappop, heappush

def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, src, dest = move
    new_stacks = [list(stack) for stack in stacks]  # Deep copy
    new_stacks[src].pop()
    new_stacks[dest].append(block)
    return tuple(tuple(stack) for stack in new_stacks)

def heuristic(stacks, goal):
    # Improved heuristic: sum of distances of each block from its goal position
    distance = 0
    goal_positions = {block: (i, j) for i, stack in enumerate(goal) for j, block in enumerate(stack)}
    for i, stack in enumerate(stacks):
        for j, block in enumerate(stack):
            if block in goal_positions:
                goal_i, goal_j = goal_positions[block]
                distance += abs(i - goal_i) + abs(j - goal_j)
    return distance

def a_star(initial_state, goal_state):
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

        for move in get_possible_moves(current_state):
            new_state = apply_move(current_state, move)
            if new_state not in visited:
                new_path = path + [move]
                cost = len(new_path) + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, new_path))

    return None

# Test with a smaller problem first
initial_state = (('A', 'G'), ('F', 'E'), ('C',), ('H', 'J', 'I', 'B', 'D'))
goal_state = (('A', 'B'), ('C', 'E', 'G', 'J'), ('D', 'F', 'H', 'I'))

solution = a_star(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(f"Move {block} from {src+1} to {dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")