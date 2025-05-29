from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Calculate a more informed heuristic
    misplaced = 0
    for stack, goal_stack in zip(state, goal_state):
        for block, goal_block in zip(stack, goal_stack):
            if block != goal_block:
                misplaced += 1
    return misplaced

def generate_moves(state):
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

def blocksworld_astar(initial_state, goal_state):
    # Priority queue for A* search
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = {initial_state: 0}

    while open_set:
        _, current_state, path = heappop(open_set)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Generate all possible moves from the current state
        for new_state, move in generate_moves(current_state):
            new_path = path + [move]
            new_cost = len(new_path)
            if new_state not in visited or new_cost < visited[new_state]:
                visited[new_state] = new_cost
                cost = new_cost + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, new_path))

# Define the initial and goal states
initial_state = (('J',), ('C', 'K', 'F', 'E', 'B', 'I'), ('G', 'A'), ('H', 'D'))
goal_state = (('F', 'G', 'H', 'I', 'J'), ('C', 'D', 'E', 'K'), ('A', 'B'))

# Find the sequence of moves
moves = blocksworld_astar(initial_state, goal_state)

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")